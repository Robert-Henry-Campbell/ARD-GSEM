# Models and summaries for the PRS evaluation. Every effect is per SD of the (within-group
# standardised) score. Count outcomes: negative binomial (MASS::glm.nb), IRR per SD, likelihood-
# ratio test against the covariates-only model, Cox-Snell pseudo-R2 of that comparison; quasi-
# Poisson fallback when theta does not converge. Mortality: Cox with age as the timescale.

.pe_na_row <- function(model, n = NA_integer_, n_events = NA_integer_, note = "") {
  data.table(model = model, n = n, n_events = n_events, estimate = NA_real_, se = NA_real_, p = NA_real_,
             effect = NA_real_, ci_low = NA_real_, ci_high = NA_real_, lr_chisq = NA_real_, lr_p = NA_real_,
             pseudo_r2 = NA_real_, theta = NA_real_, ph_p = NA_real_, note = note)
}

# y: non-negative integer counts; prs: standardised score; X: data.frame of covariates.
fit_count <- function(y, prs, X, min_events = 20L) {
  ok <- stats::complete.cases(y, prs, X)
  y <- y[ok]; prs <- prs[ok]; X <- X[ok, , drop = FALSE]
  n <- length(y); n_events <- sum(y > 0)
  if (n < 50L || n_events < min_events) return(.pe_na_row("nb", n, n_events, "too few participants or events"))
  d1 <- data.frame(y = y, prs = prs, X); d0 <- data.frame(y = y, X)
  notes <- character(0)
  fit_nb <- function(d) withCallingHandlers(
    tryCatch(MASS::glm.nb(y ~ ., data = d, control = stats::glm.control(maxit = 100)), error = function(e) e),
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") })
  m1 <- fit_nb(d1); m0 <- fit_nb(d0)
  if (inherits(m1, "error") || inherits(m0, "error")) {
    m1 <- stats::glm(y ~ ., data = d1, family = stats::quasipoisson())
    est <- stats::coef(m1)[["prs"]]; se <- sqrt(stats::vcov(m1)["prs", "prs"])
    return(data.table(model = "quasipoisson", n = n, n_events = n_events, estimate = est, se = se,
                      p = 2 * stats::pnorm(-abs(est / se)), effect = exp(est), ci_low = exp(est - 1.96 * se),
                      ci_high = exp(est + 1.96 * se), lr_chisq = NA_real_, lr_p = NA_real_, pseudo_r2 = NA_real_,
                      theta = NA_real_, ph_p = NA_real_, note = "negative binomial did not converge; quasi-Poisson (Wald p)"))
  }
  est <- stats::coef(m1)[["prs"]]; se <- sqrt(stats::vcov(m1)["prs", "prs"])
  lr <- as.numeric(2 * (stats::logLik(m1) - stats::logLik(m0)))
  data.table(model = "nb", n = n, n_events = n_events, estimate = est, se = se, p = 2 * stats::pnorm(-abs(est / se)),
             effect = exp(est), ci_low = exp(est - 1.96 * se), ci_high = exp(est + 1.96 * se),
             lr_chisq = lr, lr_p = stats::pchisq(lr, 1, lower.tail = FALSE), pseudo_r2 = 1 - exp(-lr / n),
             theta = m1$theta, ph_p = NA_real_, note = paste(unique(notes), collapse = "; "))
}

# Cox model with age as the timescale: entry = age at recruitment, exit = age at event/censoring.
fit_cox_age <- function(entry, exit, event, prs, X, min_events = 20L) {
  if (!requireNamespace("survival", quietly = TRUE)) return(.pe_na_row("cox", note = "survival not installed"))
  ok <- stats::complete.cases(entry, exit, event, prs, X) & exit > entry
  d <- data.frame(entry = entry[ok], exit = exit[ok], event = event[ok], prs = prs[ok], X[ok, , drop = FALSE])
  n <- nrow(d); n_events <- sum(d$event == 1)
  if (n < 50L || n_events < min_events) return(.pe_na_row("cox", n, n_events, "too few participants or events"))
  f <- stats::as.formula(paste("survival::Surv(entry, exit, event) ~ prs +", paste(names(X), collapse = " + ")))
  m1 <- survival::coxph(f, data = d)
  f0 <- stats::as.formula(paste("survival::Surv(entry, exit, event) ~", paste(names(X), collapse = " + ")))
  m0 <- survival::coxph(f0, data = d)
  est <- stats::coef(m1)[["prs"]]; se <- sqrt(stats::vcov(m1)["prs", "prs"])
  lr <- as.numeric(2 * (m1$loglik[2] - m0$loglik[2]))
  ph <- tryCatch(survival::cox.zph(m1)$table["prs", "p"], error = function(e) NA_real_)
  data.table(model = "cox", n = n, n_events = n_events, estimate = est, se = se, p = 2 * stats::pnorm(-abs(est / se)),
             effect = exp(est), ci_low = exp(est - 1.96 * se), ci_high = exp(est + 1.96 * se),
             lr_chisq = lr, lr_p = stats::pchisq(lr, 1, lower.tail = FALSE), pseudo_r2 = 1 - exp(-lr / n),
             theta = NA_real_, ph_p = ph, note = if (is.finite(ph) && ph < 0.05) "proportional-hazards test p < 0.05" else "")
}

# Mean outcome by score decile with normal-approximation 95% CIs (descriptive, for figures).
# Deciles are rank-based (ties broken by order): a score with few distinct values (a C+T set with
# a handful of SNPs) has coinciding quantile breaks, on which cut() fails.
decile_table <- function(prs, y) {
  ok <- is.finite(prs) & is.finite(y)
  p <- prs[ok]
  if (!length(p)) return(data.table(decile = integer(), n = integer(), mean = numeric(), se = numeric(), ci_low = numeric(), ci_high = numeric()))
  dec <- as.integer(ceiling(10 * rank(p, ties.method = "first") / length(p)))
  data.table(decile = dec, y = y[ok])[, .(n = .N, mean = mean(y), se = stats::sd(y) / sqrt(.N)), by = decile][order(decile)][
    , `:=`(ci_low = mean - 1.96 * se, ci_high = mean + 1.96 * se)][]
}

# C+T threshold per (weight-set name, group) = the threshold with the largest LR chi-square in the
# OTHER group (no tuning on the evaluation sample). Selection is per set name, not per factor:
# comparator sets are scored on the ARD count and carry factor "ARD" too. `res` needs group, name,
# kind, factor, method, threshold, lr_chisq.
ct_select_threshold <- function(res) {
  ct <- res[method == "ct" & !is.na(lr_chisq)]
  groups <- unique(ct$group)
  rbindlist(lapply(groups, function(g) {
    other <- ct[group != g]
    if (!nrow(other)) return(NULL)
    best <- other[, .SD[which.max(lr_chisq)], by = name][, .(name, kind, factor, threshold, selected_in = group, lr_chisq_other = lr_chisq)]
    best[, group := g]
    best[]
  }))
}

# Inverse-variance meta-analysis of two estimates (log scale) with Cochran's Q.
ivw_meta <- function(est1, se1, est2, se2) {
  w1 <- 1 / se1^2; w2 <- 1 / se2^2
  est <- (est1 * w1 + est2 * w2) / (w1 + w2); se <- sqrt(1 / (w1 + w2))
  q <- w1 * (est1 - est)^2 + w2 * (est2 - est)^2
  data.table(estimate = est, se = se, p = 2 * stats::pnorm(-abs(est / se)), q = q, q_p = stats::pchisq(q, 1, lower.tail = FALSE))
}
