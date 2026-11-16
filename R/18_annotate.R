# Annotate stage: variant and gene annotation for every full-QC SNP inside a locus (so every
# future credible-set member is covered) plus all locus leads.
#  - GENCODE v19 protein-coding genes: nearest TSS, overlapping gene body, and every gene
#    within candidate_window_kb of each lead (the locus candidate list)
#  - Pan-UKB nearest_genes and the Neale VEP consequence (precomputed for UKB variants)
#  - Ensembl VEP (GuideTutorial "Annotation by ANNOVAR & VEP"): offline GRCh37 cache in an
#    apptainer image, when installed; ANNOVAR (refGene, ...) when its licensed tarball is
#    present. Both are skipped with a log line otherwise.

# Protein-coding genes from a GENCODE GTF: gene_id, symbol, CHR (X = 23), start, end, strand, TSS.
read_gencode_genes <- function(gtf) {
  g <- fread(cmd = paste("zcat", shQuote(gtf), "| awk -F'\\t' '$3==\"gene\"'"), header = FALSE,
             sep = "\t", select = c(1, 4, 5, 7, 9), col.names = c("chr", "start", "end", "strand", "attr"))
  g[, gene_type := sub('.*gene_type "([^"]+)".*', "\\1", attr)]
  g <- g[gene_type == "protein_coding"]
  g[, `:=`(gene_id = sub('.*gene_id "([^".]+)[^"]*".*', "\\1", attr),
           symbol = sub('.*gene_name "([^"]+)".*', "\\1", attr),
           CHR = suppressWarnings(as.integer(sub("^chr", "", fifelse(chr == "chrX", "chr23", chr)))))]
  g <- g[!is.na(CHR)]
  g[, TSS := fifelse(strand == "+", start, end)]
  g[, .(gene_id, symbol, CHR, start, end, strand, TSS)]
}

# Nearest TSS (rolling join) and overlapping genes for variants (CHR, BP).
nearest_genes_gencode <- function(v, genes) {
  tss <- genes[, .(CHR, TSS, symbol, gene_id)]
  setkey(tss, CHR, TSS)
  q <- v[, .(row = .I, CHR, TSS = BP)]
  nn <- tss[q, on = .(CHR, TSS), roll = "nearest"]
  near <- data.table(row = nn$row, nearest_tss_gene = nn$symbol, nearest_tss_id = nn$gene_id)
  near[, nearest_tss_dist := abs(genes$TSS[match(nearest_tss_id, genes$gene_id)] - v$BP[row])]
  ov <- genes[v[, .(row = .I, CHR, BP)], on = .(CHR, start <= BP, end >= BP), nomatch = NULL,
              .(row = i.row, symbol)][, .(overlapping_genes = paste(unique(symbol), collapse = ",")), by = row]
  out <- merge(near, ov, by = "row", all.x = TRUE)
  out[order(row)]
}

# Ensembl VEP via apptainer on a VCF of the variant set; returns the tab output parsed.
run_vep <- function(config, v, work, stage) {
  sif <- config$tools$vep_sif %||% ""
  cache <- config$refs$vep_cache %||% ""
  if (!nzchar(sif) || !file.exists(sif) || !nzchar(cache) || !dir.exists(cache)) {
    log_info(stage, "VEP not installed yet (tools$vep_sif / refs$vep_cache); skipping VEP")
    return(NULL)
  }
  cfg <- config$annotate %||% list()
  vcf <- file.path(work, "vep_input.vcf")
  writeLines("##fileformat=VCFv4.2", vcf)
  fwrite(v[, .(`#CHROM` = fifelse(CHR == 23L, "X", as.character(CHR)), POS = BP, ID = SNP, REF, ALT,
               QUAL = ".", FILTER = ".", INFO = ".")][order(`#CHROM`, POS)],
         vcf, sep = "\t", append = TRUE, col.names = TRUE)
  out <- file.path(work, "vep.tsv")
  args <- c("exec", "-B", "/mnt/sdg", sif, "vep", "--input_file", vcf, "--output_file", out,
            "--tab", "--offline", "--cache", "--dir_cache", cache, "--assembly", "GRCh37",
            "--cache_version", as.character(cfg$vep_cache_version %||% 115), "--force_overwrite",
            "--fork", as.character(cfg$vep_fork %||% 16), "--buffer_size", "20000",
            "--symbol", "--canonical", "--biotype", "--numbers", "--sift", "b", "--polyphen", "b",
            "--af_gnomade", "--af_gnomadg", "--nearest", "symbol", "--distance", "5000",
            "--pick_allele_gene", "--no_stats")
  run_tool("apptainer", args, stage, file.path(work, "vep.log"), expect = out)
  lines <- readLines(out)
  hdr <- tail(grep("^#", lines, value = TRUE), 1)
  fread(text = c(sub("^#", "", hdr), grep("^#", lines, value = TRUE, invert = TRUE)), sep = "\t",
        na.strings = c("-", ""))
}

run_annotate <- function(config, sex) {
  stage <- "annotate"
  log_info(stage, sprintf("=== Annotation: %s ===", sex))
  cfg <- config$annotate %||% list()
  out_dir <- stage_dir(config, sex, stage); work <- file.path(out_dir, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  win <- as.integer(cfg$candidate_window_kb %||% 500L) * 1000L
  regions <- fread(file.path(stage_dir(config, sex, "loci"), "regions.tsv"))
  annot <- readRDS(file.path(stage_dir(config, sex, "qc_filters"), "variant_annot.rds"))
  setDT(annot)

  # variant set: full-QC SNPs inside any locus, any factor
  sets <- list()
  for (f in postgwas_factors(config)) {
    r <- regions[factor == f]
    if (!nrow(r)) next
    d <- load_postgwas(config, sex, f, "full", cols = c("SNP", "CHR", "BP"))
    hit <- d[r[, .(CHR = chr, start, end)], on = .(CHR, BP >= start, BP <= end), nomatch = NULL, .(SNP = x.SNP)]
    sets[[f]] <- hit
  }
  snps <- unique(rbindlist(sets)$SNP)
  v <- annot[SNP %in% snps]
  v[, `:=`(REF = fifelse(a1_is_ref, A1, A2), ALT = fifelse(a1_is_ref, A2, A1))]
  v <- v[!is.na(REF)]
  log_info(stage, sprintf("%s variants in %d loci", format(nrow(v), big.mark = ","), nrow(regions)))

  genes <- read_gencode_genes(config$refs$gencode_gtf)
  ng <- nearest_genes_gencode(v, genes)
  v <- cbind(v, ng[, !"row"])

  # Neale UKB VEP consequence (precomputed), matched on chr/pos/ref/alt
  nl <- fread(cmd = paste("zcat", shQuote(config$refs$neale_variants)),
              select = c("chr", "pos", "ref", "alt", "consequence", "consequence_category"),
              colClasses = list(character = "chr"), showProgress = FALSE)
  nl[, CHR := suppressWarnings(as.integer(fifelse(chr == "X", "23", chr)))]
  v <- merge(v, nl[, .(CHR, BP = pos, REF = ref, ALT = alt, neale_consequence = consequence,
                      neale_consequence_category = consequence_category)],
             by = c("CHR", "BP", "REF", "ALT"), all.x = TRUE, sort = FALSE)
  rm(nl); gc(verbose = FALSE)

  vep <- run_vep(config, v, work, stage)
  if (!is.null(vep)) {
    keep <- intersect(c("Uploaded_variation", "Consequence", "IMPACT", "SYMBOL", "Gene", "BIOTYPE",
                        "CANONICAL", "SIFT", "PolyPhen", "gnomADe_AF", "gnomADg_AF", "NEAREST", "DISTANCE"),
                      names(vep))
    vep <- vep[, ..keep]
    setnames(vep, setdiff(keep, "Uploaded_variation"), paste0("vep_", setdiff(keep, "Uploaded_variation")))
    v <- merge(v, vep[!duplicated(Uploaded_variation)], by.x = "SNP", by.y = "Uploaded_variation",
               all.x = TRUE, sort = FALSE)
  }
  if (isTRUE(cfg$annovar %||% FALSE)) {
    log_warn(stage, "annotate$annovar is TRUE but ANNOVAR support runs only once tools/annovar exists")
  }
  setorder(v, CHR, BP)
  fwrite(v, file.path(out_dir, "variants_annotated.tsv.gz"), sep = "\t")

  # per-locus leads: nearest genes and every protein-coding gene within the window
  leads <- regions[, .(region_id, factor, CHR = chr, BP = lead_bp, SNP = lead_snp, heterogeneous, mhc)]
  lead_ann <- merge(leads, v[, .(SNP, nearest_tss_gene, nearest_tss_dist, overlapping_genes, nearest_genes,
                                 neale_consequence)], by = "SNP", all.x = TRUE)
  fwrite(lead_ann[order(factor, CHR, BP)], file.path(out_dir, "lead_nearest_genes.tsv"), sep = "\t")
  cand <- genes[leads[, .(region_id, factor, lead_snp = SNP, CHR, lo = BP - win, hi = BP + win, BP)],
                on = .(CHR, TSS >= lo, TSS <= hi), nomatch = NULL,
                .(region_id, factor, lead_snp, symbol, gene_id, CHR, TSS = x.TSS, dist_to_lead = abs(x.TSS - BP))]
  fwrite(cand[order(region_id, dist_to_lead)], file.path(out_dir, "locus_candidate_genes.tsv"), sep = "\t")
  log_info(stage, sprintf("annotated %s variants; %d candidate genes across %d loci (VEP %s)",
                          format(nrow(v), big.mark = ","), nrow(cand), nrow(leads),
                          if (is.null(vep)) "pending" else "done"))
  clear_postgwas_cache()
  write_stage_manifest(stage, sex, config, list.files(out_dir, full.names = TRUE))
  invisible(v)
}
