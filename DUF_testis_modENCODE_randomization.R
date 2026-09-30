# ============================================================
# DUF TESTIS ENRICHMENT USING modENCODE RNA-seq
# Drosophila melanogaster
#
# Purpose:
#   Test whether DUF-containing genes are enriched for
#   testis-biased expression compared with random sets of
#   protein-coding Drosophila genes.
#
# Inputs:
#   1. DUF_Screen_updated_092026.xlsx
#   2. gene_rpkm_report_fb_2026_03.tsv.gz
#   3. fbgn_fbtr_fbpp_fb_2026_03.tsv.gz
#
# Main analysis:
#   - 477 DUF-containing genes
#   - modENCODE tissue RNA-seq
#   - 29 tissue samples
#   - adult male testis sample:
#       mE_mRNA_A_MateM_4d_testis
#   - 10,000 random protein-coding gene sets
# ============================================================

packages <- c("readxl","data.table","dplyr","stringr","R.utils","ggplot2")
new_packages <- packages[!(packages %in% rownames(installed.packages()))]
if (length(new_packages) > 0) install.packages(new_packages)

library(readxl)
library(data.table)
library(dplyr)
library(stringr)
library(R.utils)
library(ggplot2)

# -----------------------------
# User settings
# -----------------------------
INPUT_DIR <- "~/R"

DUF_FILE <- file.path(INPUT_DIR, "DUF_Screen_updated_092026.xlsx")
RPKM_GZ <- file.path(INPUT_DIR, "gene_rpkm_report_fb_2026_03.tsv.gz")
RPKM_TSV <- file.path(INPUT_DIR, "gene_rpkm_report_fb_2026_03.tsv")
MAP_GZ <- file.path(INPUT_DIR, "fbgn_fbtr_fbpp_fb_2026_03.tsv.gz")
MAP_TSV <- file.path(INPUT_DIR, "fbgn_fbtr_fbpp_fb_2026_03.tsv")
RESULTS_DIR <- file.path(INPUT_DIR, "DUF_testis_results")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

TESTIS_SAMPLE <- "mE_mRNA_A_MateM_4d_testis"
N_RANDOM <- 10000
RANDOM_SEED <- 12345

# -----------------------------
# Check files
# -----------------------------
cat("\nChecking input files...\n")
cat("DUF workbook:", file.exists(DUF_FILE), "\n")
cat("RNA-seq .gz:", file.exists(RPKM_GZ), "\n")
cat("FBgn-FBtr-FBpp .gz:", file.exists(MAP_GZ), "\n")

if (!file.exists(DUF_FILE)) stop("DUF workbook was not found.")
if (!file.exists(RPKM_GZ) && !file.exists(RPKM_TSV)) stop("FlyBase RNA-seq file was not found.")
if (!file.exists(MAP_GZ) && !file.exists(MAP_TSV)) stop("FlyBase FBgn-FBtr-FBpp mapping file was not found.")

# -----------------------------
# Decompress files if needed
# -----------------------------
if (!file.exists(RPKM_TSV)) {
  R.utils::gunzip(RPKM_GZ, destname = RPKM_TSV, remove = FALSE, overwrite = TRUE)
}

if (!file.exists(MAP_TSV)) {
  R.utils::gunzip(MAP_GZ, destname = MAP_TSV, remove = FALSE, overwrite = TRUE)
}

# -----------------------------
# Read FlyBase RNA-seq
# -----------------------------
cat("\nReading FlyBase RNA-seq data...\n")
rna <- data.table::fread(RPKM_TSV, skip = 5, header = TRUE)

cat("RNA-seq dimensions:", nrow(rna), "rows x", ncol(rna), "columns\n")

# -----------------------------
# Keep modENCODE tissue RNA-seq
# -----------------------------
modencode <- rna %>%
  dplyr::filter(`Parent_library_FBlc#` == "FBlc0000206")

cat("modENCODE genes:", dplyr::n_distinct(modencode$`FBgn#`), "\n")
cat("modENCODE tissues:", dplyr::n_distinct(modencode$RNASource_name), "\n")

if (nrow(modencode %>% dplyr::filter(RNASource_name == TESTIS_SAMPLE)) == 0) {
  stop("The expected modENCODE testis sample was not found.")
}

# -----------------------------
# Calculate expression metrics
# -----------------------------
testis_values <- modencode %>%
  dplyr::filter(RNASource_name == TESTIS_SAMPLE) %>%
  dplyr::transmute(
    FBgn = `FBgn#`,
    GeneSymbol = GeneSymbol,
    Testis_RPKM = as.numeric(RPKM_value)
  )

other_tissues <- modencode %>%
  dplyr::filter(RNASource_name != TESTIS_SAMPLE) %>%
  dplyr::transmute(
    FBgn = `FBgn#`,
    GeneSymbol = GeneSymbol,
    RPKM_value = as.numeric(RPKM_value)
  ) %>%
  dplyr::group_by(FBgn, GeneSymbol) %>%
  dplyr::summarise(
    Median_other_RPKM = median(RPKM_value, na.rm = TRUE),
    .groups = "drop"
  )

top_tissue <- modencode %>%
  dplyr::transmute(
    FBgn = `FBgn#`,
    GeneSymbol = GeneSymbol,
    RNASource_name = RNASource_name,
    RPKM_value = as.numeric(RPKM_value)
  ) %>%
  dplyr::group_by(FBgn, GeneSymbol) %>%
  dplyr::slice_max(order_by = RPKM_value, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    FBgn,
    GeneSymbol,
    Highest_tissue = RNASource_name,
    Highest_RPKM = RPKM_value
  )

gene_expression <- testis_values %>%
  dplyr::left_join(other_tissues, by = c("FBgn", "GeneSymbol")) %>%
  dplyr::left_join(top_tissue, by = c("FBgn", "GeneSymbol")) %>%
  dplyr::mutate(
    Testis_enrichment = log2((Testis_RPKM + 1) / (Median_other_RPKM + 1)),
    Testis_is_highest = Highest_tissue == TESTIS_SAMPLE,
    Testis_2x_enriched = Testis_enrichment >= 1 & Testis_RPKM > 0
  )

# -----------------------------
# Read DUF gene list
# -----------------------------
duf_raw <- readxl::read_excel(
  DUF_FILE,
  sheet = "DUF_genes_AF_IDs_partial_list_0"
)

duf_unique <- duf_raw %>%
  dplyr::select(Gene_name) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    Original_Gene = Gene_name,
    Clean_Gene = stringr::str_remove(Gene_name, "-R[A-Z]+$"),
    Embedded_FBgn = stringr::str_extract(Gene_name, "FBgn[0-9]+"),
    CG_Gene = stringr::str_extract(Gene_name, "CG[0-9]+")
  )

cat("Unique DUF genes:", nrow(duf_unique), "\n")

# -----------------------------
# Read FlyBase ID validation sheet
# -----------------------------
id_lookup <- readxl::read_excel(
  DUF_FILE,
  sheet = "DUF_Xmax expression"
)

id_lookup_clean <- id_lookup %>%
  dplyr::select(
    submitted_item = `#submitted_item`,
    validated_id,
    current_symbol
  ) %>%
  dplyr::filter(grepl("^FBgn", validated_id)) %>%
  dplyr::mutate(
    submitted_clean = stringr::str_remove(submitted_item, "-R[A-Z]+$")
  ) %>%
  dplyr::distinct(submitted_clean, .keep_all = TRUE)

# -----------------------------
# Match DUF genes to FlyBase IDs
# -----------------------------
expression_lookup <- gene_expression %>%
  dplyr::select(FBgn, GeneSymbol) %>%
  dplyr::distinct()

duf_match <- duf_unique %>%
  dplyr::left_join(
    id_lookup_clean %>%
      dplyr::select(submitted_clean, validated_id, current_symbol),
    by = c("Clean_Gene" = "submitted_clean")
  )

direct_lookup <- expression_lookup %>%
  dplyr::transmute(
    Clean_Gene = GeneSymbol,
    Direct_FBgn = FBgn
  )

duf_match <- duf_match %>%
  dplyr::left_join(direct_lookup, by = "Clean_Gene") %>%
  dplyr::mutate(
    FBgn_final = dplyr::coalesce(Embedded_FBgn, validated_id, Direct_FBgn)
  )

cat(
  "DUF genes with FlyBase ID:",
  sum(!is.na(duf_match$FBgn_final)),
  "of",
  nrow(duf_match),
  "\n"
)

# -----------------------------
# Attach modENCODE expression
# -----------------------------
duf_testis <- duf_match %>%
  dplyr::left_join(
    gene_expression,
    by = c("FBgn_final" = "FBgn")
  )

cat(
  "DUF genes with modENCODE expression:",
  sum(!is.na(duf_testis$Testis_RPKM)),
  "of",
  nrow(duf_testis),
  "\n"
)

# -----------------------------
# Observed DUF statistics
# -----------------------------
observed_median <- median(duf_testis$Testis_enrichment, na.rm = TRUE)
observed_mean_enrichment <- mean(duf_testis$Testis_enrichment, na.rm = TRUE)
observed_highest_n <- sum(duf_testis$Testis_is_highest, na.rm = TRUE)
observed_highest_pct <- mean(duf_testis$Testis_is_highest, na.rm = TRUE) * 100
observed_2x_n <- sum(duf_testis$Testis_2x_enriched, na.rm = TRUE)
observed_2x_pct <- mean(duf_testis$Testis_2x_enriched, na.rm = TRUE) * 100

cat("\nObserved DUF median testis enrichment:", round(observed_median, 3), "\n")
cat("Observed DUF mean testis enrichment:", round(observed_mean_enrichment, 3), "\n")
cat("Testis = highest tissue:", observed_highest_n, "\n")
cat("Percent testis-highest:", round(observed_highest_pct, 2), "%\n")
cat(">=2-fold testis enriched:", observed_2x_n, "\n")
cat("Percent >=2-fold enriched:", round(observed_2x_pct, 2), "%\n")

# -----------------------------
# Read FBgn-FBtr-FBpp mapping
# -----------------------------
fb_map <- data.table::fread(
  MAP_TSV,
  skip = "FBgn",
  header = TRUE
)

fbgn_col <- "## FlyBase_FBgn"
fbpp_col <- "FlyBase_FBpp"

protein_coding_fbgns <- fb_map %>%
  dplyr::filter(grepl("^FBpp", .data[[fbpp_col]])) %>%
  dplyr::pull(dplyr::all_of(fbgn_col)) %>%
  unique()

cat("\nProtein-coding genes identified:", length(protein_coding_fbgns), "\n")
cat(
  "DUF genes confirmed protein-coding:",
  sum(duf_testis$FBgn_final %in% protein_coding_fbgns),
  "of",
  nrow(duf_testis),
  "\n"
)

# -----------------------------
# Protein-coding background
# -----------------------------
duf_fbgns <- unique(
  duf_testis$FBgn_final[
    !is.na(duf_testis$Testis_RPKM)
  ]
)

background_pc <- gene_expression %>%
  dplyr::filter(
    FBgn %in% protein_coding_fbgns,
    !FBgn %in% duf_fbgns,
    !is.na(Testis_enrichment),
    !is.na(Testis_is_highest),
    !is.na(Testis_2x_enriched)
  )

cat("Protein-coding randomization background:", nrow(background_pc), "genes\n")

# -----------------------------
# 10,000 randomizations
# -----------------------------
set.seed(RANDOM_SEED)

n_duf <- sum(!is.na(duf_testis$Testis_RPKM))

random_results_pc <- data.frame(
  iteration = 1:N_RANDOM,
  Testis_highest_pct = NA_real_,
  Mean_enrichment = NA_real_,
  Testis_2x_pct = NA_real_
)

cat("\nRunning", N_RANDOM, "protein-coding randomizations...\n")

for (i in 1:N_RANDOM) {

  idx <- sample(
    seq_len(nrow(background_pc)),
    size = n_duf,
    replace = FALSE
  )

  random_sample <- background_pc[idx, ]

  random_results_pc$Testis_highest_pct[i] <-
    mean(random_sample$Testis_is_highest, na.rm = TRUE) * 100

  random_results_pc$Mean_enrichment[i] <-
    mean(random_sample$Testis_enrichment, na.rm = TRUE)

  random_results_pc$Testis_2x_pct[i] <-
    mean(random_sample$Testis_2x_enriched, na.rm = TRUE) * 100
}

cat("Randomization finished.\n")

# -----------------------------
# Empirical p-values
# -----------------------------
p_highest_pc <- (
  sum(random_results_pc$Testis_highest_pct >= observed_highest_pct) + 1
) / (N_RANDOM + 1)

p_mean_pc <- (
  sum(random_results_pc$Mean_enrichment >= observed_mean_enrichment) + 1
) / (N_RANDOM + 1)

p_2x_pc <- (
  sum(random_results_pc$Testis_2x_pct >= observed_2x_pct) + 1
) / (N_RANDOM + 1)

random_highest_mean_pc <- mean(random_results_pc$Testis_highest_pct)
random_enrichment_mean_pc <- mean(random_results_pc$Mean_enrichment)
random_2x_mean_pc <- mean(random_results_pc$Testis_2x_pct)
max_random_highest <- max(random_results_pc$Testis_highest_pct)

# -----------------------------
# Print final results
# -----------------------------
cat("\n=========================================\n")
cat("FINAL PROTEIN-CODING RANDOMIZATION RESULTS\n")
cat("=========================================\n")

cat("Protein-coding genes identified:", length(protein_coding_fbgns), "\n")
cat("DUF genes confirmed protein-coding:",
    sum(duf_testis$FBgn_final %in% protein_coding_fbgns), "\n")
cat("Protein-coding randomization background:", nrow(background_pc), "\n\n")

cat("Observed DUF testis-highest:", round(observed_highest_pct, 2), "%\n")
cat("Random expectation:", round(random_highest_mean_pc, 2), "%\n")
cat("Maximum random value:", round(max_random_highest, 2), "%\n")
cat("Empirical p-value:", signif(p_highest_pc, 4), "\n\n")

cat("Observed DUF mean enrichment:", round(observed_mean_enrichment, 3), "\n")
cat("Random expectation:", round(random_enrichment_mean_pc, 3), "\n")
cat("Empirical p-value:", signif(p_mean_pc, 4), "\n\n")

cat("Observed DUF >=2-fold:", round(observed_2x_pct, 2), "%\n")
cat("Random expectation:", round(random_2x_mean_pc, 2), "%\n")
cat("Empirical p-value:", signif(p_2x_pc, 4), "\n")

# -----------------------------
# Save result tables
# -----------------------------
write.csv(
  duf_testis,
  file.path(RESULTS_DIR, "DUF_gene_testis_expression_results.csv"),
  row.names = FALSE
)

write.csv(
  random_results_pc,
  file.path(RESULTS_DIR, "DUF_testis_randomization_10000_protein_coding.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  Metric = c(
    "DUF genes analyzed",
    "Protein-coding genes identified",
    "Protein-coding background genes",
    "Observed testis-highest percent",
    "Random testis-highest mean percent",
    "Maximum random testis-highest percent",
    "Testis-highest empirical p-value",
    "Observed mean enrichment",
    "Random mean enrichment",
    "Mean-enrichment empirical p-value",
    "Observed >=2-fold percent",
    "Random >=2-fold mean percent",
    ">=2-fold empirical p-value"
  ),
  Value = c(
    n_duf,
    length(protein_coding_fbgns),
    nrow(background_pc),
    observed_highest_pct,
    random_highest_mean_pc,
    max_random_highest,
    p_highest_pc,
    observed_mean_enrichment,
    random_enrichment_mean_pc,
    p_mean_pc,
    observed_2x_pct,
    random_2x_mean_pc,
    p_2x_pc
  )
)

write.csv(
  summary_table,
  file.path(RESULTS_DIR, "DUF_testis_randomization_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# Final figure
# -----------------------------
p_testis_pc <- ggplot(
  random_results_pc,
  aes(x = Testis_highest_pct)
) +
  geom_histogram(
    bins = 40,
    fill = "grey75",
    color = "black",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = observed_highest_pct,
    color = "blue",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = observed_highest_pct,
    y = Inf,
    label = paste0("DUF genes = ", round(observed_highest_pct, 1), "%"),
    hjust = 1.05,
    vjust = 2,
    size = 4
  ) +
  annotate(
    "text",
    x = observed_highest_pct,
    y = Inf,
    label = paste0(
      "Empirical p = ",
      format(p_highest_pc, scientific = TRUE, digits = 2)
    ),
    hjust = 1.05,
    vjust = 4,
    size = 4
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "DUF-Containing Genes Are Enriched for Testis Expression",
    subtitle = paste0(
      format(N_RANDOM, big.mark = ","),
      " random protein-coding gene sets (n = ",
      n_duf,
      " per set)"
    ),
    x = "Genes with testis as highest-expression tissue (%)",
    y = "Number of random gene sets"
  )

print(p_testis_pc)

ggsave(
  filename = file.path(RESULTS_DIR, "DUF_testis_randomization_protein_coding.png"),
  plot = p_testis_pc,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(RESULTS_DIR, "DUF_testis_randomization_protein_coding.pdf"),
  plot = p_testis_pc,
  width = 7,
  height = 5
)

cat("\nAnalysis complete.\n")
cat("Results were saved to:\n", RESULTS_DIR, "\n")
