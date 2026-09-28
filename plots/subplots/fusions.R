library(tidyverse)

results_root <- "data/fusions"

# Display labels
run_levels <- c("Pre-dedup", "Post-dedup")

#-------------------------
# VISUALIZATION CONFIGS
#-------------------------
run_colors <- c("Pre-dedup"  = "#01476C",
                "Post-dedup" = "#8EBCDA")

detect_levels <- c("Detected", "Not detected")

my_theme <- theme_minimal(base_size = 9) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    panel.spacing    = unit(1.5, "lines")
  )

#-------------------------
# FUSION (JAFFAL) CONFIGURATION
#-------------------------
# JAFFAL run directory -> display label
fus_run_dirs <- c("Pre-dedup"  = "lb_bulk_h146_original",
                  "Post-dedup" = "lb_bulk_h146_nailpolish")

# H146 ONT dRNA JAFFAL calls as the fusion truth set. Same cell line, different library prep, and no deduplication applied to it
FUS_DRNA_CSV      <- file.path(results_root, "ref", "comprehensive_downsampled_data_fusion_result_include_ONT_cDNA_dedup.csv")
FUS_DRNA_CELLLINE <- "H146"
FUS_DRNA_PLATFORM <- "ONT dRNA"

# confidence tiers order
fus_class_levels <- c("HighConfidence", "LowConfidence", "PotentialTransSplicing")

# fusion threshold, applied to the CALLS. Figures 1 and 2 vary this.
fus_min_spanning_reads <- 3

# Truth-side support gate
fus_drna_min_spanning_reads <- 3

#-------------------------
# SNV (Clair3-RNA) CONFIGURATION
#-------------------------
# Single reference set: LongBench's own GRCh38 with the 115 spike-in contigs and 169 GL/KI decoys
snv_root <- file.path(results_root, "output", "SNV")

snv_run_dirs <- c("Pre-dedup"  = "original",
                  "Post-dedup" = "nailpolish")

# Clair3-RNA
VCF_NAME <- "output_enable_phasing.vcf.gz"

# filter only chr1-22, X,Y
MAIN_CTGS <- c(paste0("chr", 1:22), "chrX", "chrY")

# Sanger WES truth, written by snv/03b_h146_sanger_wes_truth_vcf_processing.Rmd.
SANGER_TRUTH_TSV <- file.path(snv_root, "ref", "H146_sanger_wes_truth_table.tsv")

# Per-site depth at the truth positions
DEPTH_TSVS_SANGER <- file.path(snv_root, "depth",
                               paste0("truth_depth_sanger_", snv_run_dirs, ".tsv"))
names(DEPTH_TSVS_SANGER) <- run_levels

# SNV threshold
MIN_DEPTH_BOTH <- 10

# Ready-made VAF strata carried in the truth TSV, used by the truth-set import table below
snv_vaf_levels <- c("subclonal (<0.25)", "het-like (0.25-0.90)", "hom-like (>=0.90)")

snv_vtype_levels <- c("SNV", "indel", "other")

#-------------------------
# FUSION HELPERS
#-------------------------

# Normalise a JAFFA fusion name to an orientation-independent key so that A::B and B::A collapse to one call
# NOTE: JAFFA uses "::" (not ":"), and gene symbols may themselves contain hyphens (TMEM41A-DT) or be raw Ensembl IDs (ENSG00000307038), so split on the literal "::" and do not touch the tokens beyond trim/upper
pair_key_from_fusion <- function(x) {
  map_chr(x, function(s) {
    g <- str_split_1(as.character(s), fixed("::"))
    g <- toupper(trimws(g))
    g <- g[!is.na(g) & g != ""]
    paste(sort(g), collapse = "::")
  })
}

#-------------------------
# SNV HELPERS - verbatim from the Sanger document
#-------------------------
# Read a Clair3-RNA VCF as a plain table. comment = "##" drops the meta lines but keeps the "#CHROM" line, which becomes the header, hence the rename. Everything is read as character and coerced explicitly, because QUAL can be "." and AD is comma-separated ("25,19")
read_clair3_vcf <- function(path) {
  vcf <- read_tsv(path, comment = "##", col_types = cols(.default = col_character()),
                  progress = FALSE)
  names(vcf)[1] <- "CHROM"

  # FORMAT is GT:GQ:DP:AD:AF on every record in these files. Assert it rather than assume: splitting the sample column positionally is only safe while that holds.
  stopifnot(n_distinct(vcf$FORMAT) == 1)
  fmt <- str_split_1(vcf$FORMAT[1], fixed(":"))

  parts <- str_split_fixed(vcf$SAMPLE, fixed(":"), length(fmt))
  colnames(parts) <- fmt

  bind_cols(vcf, as_tibble(parts))
}

# Coerce a character column to numeric, treating the literal string "NA" as missing
num <- function(x) suppressWarnings(as.numeric(na_if(as.character(x), "NA")))

# Sorted-allele genotype, so 1|0 and 0/1 compare equal
norm_gt <- function(gt) {
  map_chr(gt, function(g) {
    if (is.na(g) || g %in% c(".", "./.", ".|.")) return(".")
    paste(sort(str_split_1(str_replace_all(g, fixed("|"), "/"), fixed("/"))), collapse = "/")
  })
}

variant_type <- function(ref, alt) {
  case_when(nchar(ref) == 1 & nchar(alt) == 1 ~ "SNV",
            nchar(ref) != nchar(alt)          ~ "indel",
            TRUE                              ~ "other")
}

# Label a truth set against the calls
# Matching is exact on CHROM:POS:REF:ALT
label_against_truth <- function(truth_tbl, call_lookup, score_gt,
                                extra_cols = character()) {
  truth_cols <- c("key", "HugoSymbol", "ProteinChange", "vtype", "AF_truth", extra_cols)
  if (score_gt) truth_cols <- c(truth_cols, "GT_truth")

  tt <- truth_tbl %>%
    mutate(GT_truth = if (score_gt) GT else NA_character_, AF_truth = AF) %>%
    dplyr::select(all_of(truth_cols))

  levs <- if (score_gt) c("Match", "Discordant GT", "Not called")
          else          c("Called", "Not called")

  out <- expand_grid(key = truth_tbl$key,
                     run = factor(run_levels, levels = run_levels)) %>%
    left_join(tt, by = "key") %>%
    left_join(call_lookup, by = c("key", "run")) %>%
    mutate(
      status = if (score_gt) {
        case_when(is.na(GT)                      ~ "Not called",
                  gt_norm == norm_gt(GT_truth)   ~ "Match",
                  TRUE                           ~ "Discordant GT")
      } else {
        if_else(is.na(GT), "Not called", "Called")
      },
      status = factor(status, levels = levs)
    )

  stopifnot(nrow(out) == nrow(truth_tbl) * length(run_levels),
            !any(is.na(out$status)))
  out
}

# Per-site depth at a truth set's positions, widened to one row per truth variant.
# The NA guard is the important part. Each truth set has its OWN depth files, because the two sets sit at different positions. Point this at the wrong pair and most variants join to NA, covered_both becomes NA rather than FALSE, and a sum(..., na.rm = TRUE) downstream would report a small, entirely plausible number instead of failing. So: no na.rm anywhere below, and a hard stop on any missing row.
join_truth_depth <- function(truth_tbl, depth_tsvs) {
  depth <- depth_tsvs %>%
    map(read_tsv, col_types = cols(CHROM = col_character(), POS = col_double(),
                                   depth = col_double(), arm = col_character()),
        progress = FALSE) %>%
    bind_rows(.id = "run") %>%
    mutate(run = factor(run, levels = run_levels)) %>%
    dplyr::select(run, CHROM, POS, depth)

  out <- truth_tbl %>%
    dplyr::select(key, CHROM, POS) %>%
    left_join(depth, by = c("CHROM", "POS")) %>%
    dplyr::select(key, run, depth) %>%
    pivot_wider(names_from = run, values_from = depth)

  missing <- rowSums(is.na(dplyr::select(out, all_of(run_levels)))) > 0
  if (any(missing)) {
    stop(sum(missing), " of ", nrow(out), " truth positions have no depth row. The depth ",
         "files were built for a different truth set - re-run 05_coverage.sbatch with the ",
         "matching TRUTHSET. Proceeding would make covered_both NA rather than FALSE and ",
         "silently shrink the denominator.")
  }

  out %>%
    mutate(covered_both = `Pre-dedup` >= MIN_DEPTH_BOTH & `Post-dedup` >= MIN_DEPTH_BOTH)
}

#-------------------------
# THE FIGURE HELPERS
#-------------------------

# Per-arm count, in run_levels order
n_by_run <- function(d) {
  out <- d %>% count(run, name = "n", .drop = FALSE) %>% arrange(run)
  stopifnot(identical(out$run, factor(run_levels, levels = run_levels)))
  out$n
}

# Recovered / denominator per arm, from a truth-anchored frame carrying 'run' and a 'status' factor on detect_levels
recovery_tbl <- function(labelled, denom_by = NULL) {
  labelled %>%
    group_by(across(all_of(c("run", denom_by)))) %>%
    summarise(recovered = sum(status == "Detected"),
              denom     = n(),
              .groups   = "drop") %>%
    arrange(run) %>%
    mutate(label = sprintf("%d/%d\n(%.1f%%)", recovered, denom,
                           100 * recovered / denom))
}

# Truth-anchored fusion recovery at a spanning-read threshold
fus_recovery_at <- function(min_sr) {
  called <- fus_calls %>%
    dplyr::filter(spanning_reads >= min_sr) %>%
    transmute(run, feature = pair_key) %>%
    distinct() %>%
    mutate(called = TRUE)

  out <- expand_grid(feature = fus_truth_pairs_drna,
                     run     = factor(run_levels, levels = run_levels)) %>%
    left_join(called, by = c("run", "feature")) %>%
    mutate(status = factor(if_else(is.na(called), "Not detected", "Detected"),
                           levels = detect_levels))
  out
}

# measure x run, each measure indexed to its OWN pre-dedup value
fig_idx <- function(calls_label, calls_n, truth_label, truth_n) {
  out <- bind_rows(
    tibble(measure = calls_label,
           run     = factor(run_levels, levels = run_levels),
           n       = calls_n),
    tibble(measure = truth_label,
           run     = factor(run_levels, levels = run_levels),
           n       = truth_n)
  ) %>%
    group_by(measure) %>%
    mutate(pct = 100 * n / n[run == "Pre-dedup"]) %>%
    ungroup() %>%
    mutate(measure = fct_inorder(measure),
           label   = sprintf("%s\n%.0f%%", scales::comma(n), pct))

  # pct divides by this, and a zero pre-dedup bar would make the panel Inf rather than fail.
  stopifnot(all(out$n[out$run == "Pre-dedup"] > 0))
  out
}

fig_plot <- function(idx, title, subtitle, caption = NULL) {
  ggplot(idx, aes(x = measure, y = pct, fill = run)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.68) +
    # The pre-dedup bar IS 100% by construction; the rule makes the post-dedup bar readable against it without having to compare bar heights across groups
    geom_hline(yintercept = 100, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_text(aes(label = label), position = position_dodge(width = 0.8),
              vjust = -0.3, size = 2.5, lineheight = 0.9) +
    scale_fill_manual(values = run_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22)),
                       breaks = seq(0, 90, by = 30),
                       labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "Relative to pre-dedup",
         title = title) +
    my_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"))
}

#-------------------------
# PREPROCESSING
#-------------------------

fus_csv_paths <- file.path(results_root, "output", "JAFFAL", fus_run_dirs,
                           "jaffa_results.csv")
names(fus_csv_paths) <- run_levels

fus_rows <- fus_csv_paths %>%
  map(read_csv, col_types = cols(.default = col_character())) %>%
  bind_rows(.id = "run") %>%
  mutate(
    run            = factor(run, levels = run_levels),
    fusion         = `fusion genes`,
    spanning_reads = as.numeric(`spanning reads`),
    classification = factor(classification, levels = fus_class_levels),
    pair_key       = pair_key_from_fusion(fusion)
  )

fus_drna_all <- read_csv(FUS_DRNA_CSV, col_types = cols(.default = col_character()),
                         progress = FALSE)

fus_drna_h146 <- fus_drna_all %>%
  filter(CellLine == FUS_DRNA_CELLLINE, Platform == FUS_DRNA_PLATFORM) %>%
  mutate(spanning_reads = as.numeric(`spanning.reads`),
         classification = factor(classification, levels = fus_class_levels),
         pair_key       = pair_key_from_fusion(`fusion.genes`))

fus_drna_kept <- fus_drna_h146 %>%
  filter(spanning_reads >= fus_drna_min_spanning_reads)

# Rows collapse to fewer pairs because the same junction is reported on several contigs
fus_truth_pairs_drna_raw <- unique(fus_drna_kept$pair_key)

fus_truth_pairs_drna <- fus_truth_pairs_drna_raw

# One row per (run, pair): keep the maximum spanning-read support and the best (highest-confidence) classification seen for that pair
fus_calls <- fus_rows %>%
  group_by(run, pair_key) %>%
  summarise(
    fusion         = first(fusion),
    spanning_reads = max(spanning_reads),
    classification = factor(fus_class_levels[min(as.integer(classification))],
                            levels = fus_class_levels),
    n_rows         = n(),
    .groups        = "drop"
  ) %>%
  mutate(drna_pair_support = pair_key %in% fus_truth_pairs_drna)

stopifnot(!any(is.na(fus_calls$spanning_reads)),
          !any(is.na(fus_calls$classification)))


snv_vcf_paths <- file.path(snv_root, "clair3_rna", snv_run_dirs, VCF_NAME)
names(snv_vcf_paths) <- run_levels

snv_calls_raw <- snv_vcf_paths %>%
  map(read_clair3_vcf) %>%
  bind_rows(.id = "run") %>%
  mutate(
    run   = factor(run, levels = run_levels),
    POS   = as.numeric(POS),
    QUAL  = as.numeric(QUAL),
    DP    = as.numeric(DP),
    AF    = as.numeric(AF),
    vtype = factor(variant_type(REF, ALT), levels = snv_vtype_levels),
    key   = paste(CHROM, POS, REF, ALT, sep = ":")
  )

# filter PASS calls
snv_calls <- snv_calls_raw %>%
  filter(CHROM %in% MAIN_CTGS, FILTER == "PASS")


snv_truth <- read_tsv(SANGER_TRUTH_TSV, col_types = cols(.default = col_character()),
                      progress = FALSE) %>%
  mutate(POS   = as.numeric(POS),
         AF    = num(AF),
         DP    = num(DP),
         # Rebuild the key rather than trusting the one in the file, so a schema change upstream fails an assertion instead of silently mismatching
         key2  = paste(CHROM, POS, REF, ALT, sep = ":"),
         vtype = variant_type(REF, ALT),
         across(c(npgl, drv, exonic, dp_ok), ~ .x == "TRUE"),
         npgl_class = factor(npgl_class,
                             levels = c("germline panel hit", "no panel hit")),
         vaf_stratum = factor(vaf_stratum, levels = snv_vaf_levels))

stopifnot(identical(snv_truth$key, snv_truth$key2),
          !any(duplicated(snv_truth$key)),
          all(snv_truth$vtype == "SNV"),
          all(snv_truth$CHROM %in% MAIN_CTGS),
          !any(is.na(snv_truth$vaf_stratum)),
          !any(is.na(snv_truth$exonic)))

snv_truth_depth <- join_truth_depth(snv_truth, DEPTH_TSVS_SANGER)

# pivot_wider would inflate this on a duplicated CHROM:POS in the depth files
stopifnot(nrow(snv_truth_depth) == nrow(snv_truth))

snv_n_cov <- sum(snv_truth_depth$covered_both)

snv_call_lookup <- snv_calls %>%
  dplyr::select(run, key, QUAL, DP, AF, GT) %>%
  mutate(gt_norm = norm_gt(GT))

snv_labelled <- label_against_truth(
    snv_truth, snv_call_lookup, score_gt = FALSE, extra_cols = c("exonic")) %>%
  # Normalise onto the shared vocabulary so one recovery helper serves both assays
  # label_against_truth() is kept byte-identical to the source, so the recode is done here and stays visible
  mutate(status = factor(if_else(status == "Not called", "Not detected", "Detected"),
                         levels = detect_levels)) %>%
  left_join(snv_truth_depth %>% dplyr::select(key, covered_both), by = "key")

# The two SNV denominators. These are NOT nested: the covered set is whatever both arms sequence deeply, which includes intronic positions, while the exonic set is everything the truth file marks as exonic regardless of depth
snv_labelled_cov    <- snv_labelled %>% filter(covered_both)
snv_labelled_exonic <- snv_labelled %>% filter(exonic)

#-------------------------
# FIGURE 1
#-------------------------

stopifnot(fus_min_spanning_reads >= min(fus_calls$spanning_reads),
          fus_min_spanning_reads <= max(fus_calls$spanning_reads))

fus_thr_lab <- fus_recovery_at(fus_min_spanning_reads)
fus_thr_rec <- recovery_tbl(fus_thr_lab)

fus_thr_idx <- fig_idx(
  calls_label = sprintf("Fusion pairs called\n(>= %d spanning reads)",
                        fus_min_spanning_reads),
  calls_n     = n_by_run(fus_calls %>%
                           filter(spanning_reads >= fus_min_spanning_reads)),
  truth_label = sprintf("dRNA fusions recovered\n(of %d dRNA pairs,\n>= %d spanning reads)",
                        length(fus_truth_pairs_drna), fus_drna_min_spanning_reads),
  truth_n     = fus_thr_rec$recovered)

# plot
p_fus_thr <- fig_plot(
  fus_thr_idx,
  title = NULL)

#-------------------------
# FIGURE 2
#-------------------------

snv_thr_rec <- recovery_tbl(snv_labelled_cov)

stopifnot(all(snv_thr_rec$denom == snv_n_cov))

snv_thr_idx <- fig_idx(
  calls_label = sprintf("PASS calls\n(DP >= %d)", MIN_DEPTH_BOTH),
  calls_n     = n_by_run(snv_calls %>% filter(DP >= MIN_DEPTH_BOTH)),
  truth_label = sprintf("Sanger truth recovered\n(of %d sites covered >= %dx in both)",
                        snv_n_cov, MIN_DEPTH_BOTH),
  truth_n     = snv_thr_rec$recovered)

p_snv_thr <- fig_plot(
  snv_thr_idx,
  title    = sprintf("H146 ONT cDNA SNVs: at >= %dx coverage", MIN_DEPTH_BOTH))


#-------------------------
# FIGURE 3
#-------------------------

snv_all_rec <- recovery_tbl(snv_labelled_exonic)

stopifnot(all(snv_all_rec$denom == sum(snv_truth$exonic)))

snv_all_idx <- fig_idx(
  calls_label = "PASS calls\n(all, chr1-22, X, Y)",
  calls_n     = n_by_run(snv_calls),
  truth_label = sprintf("Sanger truth recovered\n(of %d exonic sites)",
                        snv_all_rec$denom[1]),
  truth_n     = snv_all_rec$recovered)

# plot
p_snv_all <- fig_plot(
  snv_all_idx,
  title    = "H146 ONT cDNA SNVs: every call")