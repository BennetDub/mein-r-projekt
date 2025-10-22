# HIPATHIA PIPELINE – TEMPLAT   44444444
# ------------------------------------------------------------
# Purpose: Run your own gene expression data through hipathia
# Author: <your name>
# Date: <today>
# ------------------------------------------------------------
# INPUT ASSUMPTIONS
# - Expression matrix: rows = genes (IDs), columns = samples. File: CSV/TSV.
# - Metadata (design): at least columns `sample` (must match matrix column names)
#   and `group` (e.g., Tumor/Normal).
# - Gene IDs: ENSEMBL or ENTREZ etc.; hipathia will map to the required IDs.
# ------------------------------------------------------------

# =====================
# 0) PARAMETER BLOCK
# =====================
# 👉 Adjust these to your data.

# 0.1 File paths
expr_file   <- "./rna_data_filtered_tumor_Hipathia.tsv"   # Genes x Samples; first column = gene IDs; header = sample IDs 
design_file <- "./cell_metadata_filtered_hipathia.tsv"              # Tab-delimited; columns: sample, group (and optional others) Metadata
output_dir  <- "./hipathia_report"         # Output folder for the HTML report

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# 0.2 Organism / Pathways
species     <- "hsa"                        # "hsa" (Human), "mmu" (Mouse), "rno" (Rat)
# Optional: load only certain KEGG IDs (saves time/memory). Otherwise: NULL
pathways_list <- NULL                    # e.g., c("hsa03320","hsa04014") or NULL for all

# 0.3 Groups (order defines direction: g1 vs g2)
group1 <- "Tumor"                           # g1 (interpreted as "up" in comparisons)
group2 <- "Normal"                          # g2

# 0.4 Normalization (see vignette section 3.2)
by_quantiles         <- FALSE               # TRUE = quantile normalization first
percentil_mode       <- FALSE               # TRUE = percentile scaling instead of linear scaling
truncation_percentil <- NULL                 # e.g., 0.95; otherwise keep NA

# 0.5 PCA / feature selection
use_pca             <- TRUE
max_pca_features    <- NA                   # Default: ncol(X); set if you want a fixed cap

# 0.6 Statistics / FDR
p_adj_method        <- "BH"                 # Benjamini–Hochberg
alpha               <- 0.05                 # FDR threshold

# 0.7 Visualization / local server
run_local_server    <- FALSE                # TRUE starts visualize_report(...) at the end
server_port         <- 4000

# =====================
# 1) SETUP
# =====================
needed <- c("hipathia","SummarizedExperiment","S4Vectors","limma")

for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (Sys.getenv("CONDA_PREFIX") != "") {
      stop("Package ", pkg, " fehlt. Bitte mit Conda installieren.")
    } else {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
}


set.seed(1)

# =====================
# 2) LOAD DATA (robust)
# =====================
read_table_guess <- function(path) {
  stopifnot(file.exists(path))
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","txt")) {
    utils::read.table(path, sep = "\t", header = TRUE, quote = "\"",
                      check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  } else {
    utils::read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

expr_df   <- read_table_guess(expr_file)
design_df <- read_table_guess(design_file)

stopifnot("sample" %in% colnames(design_df), "group" %in% colnames(design_df))
stopifnot(ncol(expr_df) >= 2)   # 1x IDs + >=1 Sample

# 1) IDs exakt einmal setzen
gene_ids <- as.character(expr_df[[1]])
stopifnot(!anyNA(gene_ids))
rownames(expr_df) <- gene_ids
expr_df[[1]] <- NULL

# 2) Samples abgleichen
samples <- intersect(colnames(expr_df), as.character(design_df$sample))
if (length(samples) == 0) stop("Kein Overlap zwischen Expression-Spalten und design$sample.")
if (length(samples) < ncol(expr_df)) {
  warning("Ignoriere ", ncol(expr_df)-length(samples), " Expression-Spalten ohne Design-Eintrag.")
}
expr_df   <- expr_df[, samples, drop = FALSE]
design_df <- design_df[match(samples, design_df$sample), , drop = FALSE]
stopifnot(all(design_df$sample == colnames(expr_df)))

# 3) Numerik hart erzwingen
expr_mat_num <- suppressWarnings(as.matrix(apply(expr_df, 2, as.numeric)))
if (anyNA(expr_mat_num)) {
  nbad <- sum(is.na(expr_mat_num))
  stop("Nicht-numerische Einträge in Expressionsmatrix (", nbad, " NAs nach coercion). Datei prüfen.")
}
expr_df <- as.data.frame(expr_mat_num, check.names = FALSE)

message("Data: ", nrow(expr_df), " genes x ", ncol(expr_df), " samples")

# =====================
# 3) ID TRANSLATION & NORMALIZATION
# =====================
expr_mat <- as.matrix(expr_df)
mode(expr_mat) <- "numeric"

# 3.1 Gene-ID-Übersetzung (Vignette 3.1)
tr <- hipathia::translate_data(expr_mat, species = species)

# Versionstolerant: Slot-Namen können variieren
`%||%` <- function(a, b) if (!is.null(a)) a else b
trans_data <- tr$exp %||% tr$matrix %||% tr$translated
if (is.null(trans_data)) stop("translate_data() lieferte kein exp/matrix/translated – hipathia-Version prüfen.")

# 3.2 Skalierung/Normalisierung (Vignette 3.2)
message("\n[Normalization] …")
exp_data <- hipathia::normalize_data(
  trans_data,
  by_quantiles = by_quantiles,
  percentil    = percentil_mode,
  truncation_percentil = truncation_percentil
)

# =====================
# 4) LOAD PATHWAYS
# =====================
message("
[Load pathways] …")
if (is.null(pathways_list)) {
  pathways <- load_pathways(species = species)
} else {
  pathways <- load_pathways(species = species, pathways_list = pathways_list)
}
message("Loaded pathways: ", length(get_pathways_list(pathways)))

# =====================
# 5) BUILD SUMMARIZEDEXPERIMENT (optional)
# =====================
# If you want a SE object for bookkeeping:
col_data <- S4Vectors::DataFrame(design_df[, setdiff(colnames(design_df), "sample"), drop = FALSE], row.names = design_df$sample)
se <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(raw = exp_data), colData = col_data)

# =====================
# 6) RUN HIPATHIA
# =====================
message("
[Hipathia] Computing activations (effector subpaths) …")
# decompose = FALSE for overview; TRUE only for targeted deep dives (many more features)
results <- hipathia(exp_data, pathways, decompose = FALSE)
print(results)

# =====================
# 7) EXTRACT FEATURES (SUBPATHS / FUNCTIONS)
# =====================
path_vals <- get_paths_data(results, matrix = TRUE)  # Matrix: subpaths x samples
# Optional function level (Uniprot/GO)
# uniprot_se <- quantify_terms(results, pathways, dbannot = "uniprot")
# go_se      <- quantify_terms(results, pathways, dbannot = "GO")

# =====================
# 8) STATS: WILCOXON g1 VS g2
# =====================
message("
[Wilcoxon] ", group1, " vs ", group2)
sample_group <- col_data$group
names(sample_group) <- rownames(col_data)
stopifnot(all(names(sample_group) == colnames(path_vals)))

comp <- do_wilcoxon(path_vals, sample_group, g1 = group1, g2 = group2)
# FDR correction
if (!"p.value" %in% colnames(comp)) {
  # Some hipathia versions name the p-value column differently
  pcol <- grep("p", colnames(comp), value = TRUE)[1]
  colnames(comp)[colnames(comp) == pcol] <- "p.value"
}
comp$q.value <- p.adjust(comp$p.value, method = p_adj_method)
# Ergebnisse persistieren (optional, aber sinnvoll)
utils::write.table(comp,
  file = file.path(output_dir, "subpath_differential.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE
)

sig <- comp[comp$q.value <= alpha, , drop = FALSE]
if (nrow(sig) == 0) {
  warning("Keine signifikanten Subpaths bei q<=", alpha)
} else {
  utils::write.table(sig,
    file = file.path(output_dir, "subpath_significant.tsv"),
    sep = "\t", quote = FALSE, row.names = TRUE
  )
}


# =====================
# 10) PCA (optional; visual QC)
# =====================
if (use_pca) {
  message("\n[PCA] …")

  # Features nach p-Wert sortieren (falls vorhanden)
  if (!is.null(comp) && "p.value" %in% colnames(comp)) {
    ranked_ids <- rownames(path_vals)[order(comp$p.value, decreasing = FALSE)]
    ranked_ids <- ranked_ids[ranked_ids %in% rownames(path_vals)]
    ranked_path_vals <- path_vals[ranked_ids, , drop = FALSE]
  } else {
    ranked_path_vals <- path_vals
  }

  nfeat  <- nrow(ranked_path_vals)
  nsamp  <- ncol(ranked_path_vals)
  k_cap  <- max(0, nsamp - 1)  # princomp braucht Variablen < Beobachtungen
  k_user <- if (is.na(max_pca_features)) nfeat else min(nfeat, max_pca_features)
  k      <- min(nfeat, k_cap, k_user)

  if (k < 2) {
    message("Zu wenige zulässige Features für PCA (k=", k, ", nsamp=", nsamp, "). Schritt übersprungen.")
  } else {
    X <- ranked_path_vals[seq_len(k), , drop = FALSE]
    pca_model <- hipathia::do_pca(X)
    try( hipathia::pca_plot(pca_model, sample_group, legend = TRUE), silent = TRUE )
  }
}


# =====================
# 11) NDE-LEVEL DE COLORS (limma) & COMPARISON PLOT
# =====================
message("\n[Node DE + pathway plot] …")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

colors_de <- node_color_per_de(results, pathways, sample_group, group1, group2)

plot_id <- if (!is.null(pathways_list)) pathways_list[1] else get_pathways_list(pathways)[1]
try({
  pathway_comparison_plot(comp, metaginfo = pathways, pathway = plot_id, node_colors = colors_de)
}, silent = TRUE)

# =====================
# 12) CREATE REPORT & SAVE SUMMARIES
# =====================
message("\n[Report] Writing HTML report into ", normalizePath(output_dir, mustWork = FALSE))

pathways_summary <- get_pathways_summary(comp, pathways)
pathways_summary$percent_significant_paths <- 100 * pathways_summary$num_significant_paths / pmax(1, pathways_summary$num_total_paths)
pathways_summary <- pathways_summary[order(-pathways_summary$percent_significant_paths, -pathways_summary$num_significant_paths), ]
utils::write.table(pathways_summary,
                   file = file.path(output_dir, "pathways_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

report <- create_report(comp, pathways, output_dir, node_colors = colors_de)

if (isTRUE(run_local_server)) {
  visualize_report(report, port = server_port)
}

# =====================
# 13) QUALITY NOTES / LOG
# =====================
message("\n[NOTES]")
message("- Check translate/normalize output: counts of untranslated/multihit IDs.")
message("- Many imputed missing genes in hipathia ⇒ interpret results with caution.")
message(paste0("- 'up'/'down' direction refers to g1=", group1, " vs g2=", group2, "."))
message("- For functional level: quantify_terms(dbannot='uniprot'/'GO').")
message(paste0("- Applied FDR (", p_adj_method, ") threshold: q<", alpha, "."))

message("\n[DONE]")
