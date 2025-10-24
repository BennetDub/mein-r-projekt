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
.log <- function(...) message(strftime(Sys.time(), "%H:%M:%S "), sprintf(...))

.time_it <- function(expr, label="step") {
  t0 <- Sys.time(); on.exit(.log("%s done in %.1f min", label, as.numeric(difftime(Sys.time(), t0, units="mins"))))
  force(expr)
}

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

# --- A: Diagnose nach dem Laden ---
cat("Expr-Dims: ", nrow(expr_df), " Gene x ", ncol(expr_df), " Samples\n", sep = "")
cat("Erste 5 Gene: ", paste(head(rownames(expr_df), 5), collapse=", "), "\n", sep = "")
cat("Vermutete ID-Form (Symbol!=nur Ziffern): ",
    any(!grepl("^\\d+$", rownames(expr_df))), "\n", sep = "")


# 3) Numerik hart erzwingen
expr_mat_num <- suppressWarnings(as.matrix(apply(expr_df, 2, as.numeric)))
if (anyNA(expr_mat_num)) {
  nbad <- sum(is.na(expr_mat_num))
  stop("Nicht-numerische Einträge in Expressionsmatrix (", nbad, " NAs nach coercion). Datei prüfen.")
}
expr_df <- as.data.frame(expr_mat_num, check.names = FALSE)

message("Data: ", nrow(expr_df), " genes x ", ncol(expr_df), " samples")

# ==== 3) ID TRANSLATION & NORMALIZATION (robust) 


# 3) ROBUSTES ID-MAPPING & NORMALISIERUNG
# =====================

species <- "hsa"  # Human

# 3.0: Gene-IDs evtl. von Anführungszeichen säubern
rownames(expr_df) <- gsub("^['\"]+|['\"]+$", "", rownames(expr_df))

# Numerik erzwingen
expr_mat0 <- as.matrix(expr_df)
storage.mode(expr_mat0) <- "double"

# Helper
`%||%` <- function(a,b) if (!is.null(a)) a else b
grab_mat <- function(tr) {
  if (is.matrix(tr)) return(tr)
  if (is.list(tr))   return(tr$exp %||% tr$matrix %||% tr$translated)
  if (inherits(tr,"SummarizedExperiment")) return(SummarizedExperiment::assay(tr,1))
  NULL
}

# 3.1: ID-Typen bestimmen
ids <- rownames(expr_mat0)
looks_numeric <- all(grepl("^[0-9]+$", head(ids, 100), perl = TRUE)) || mean(grepl("^[0-9]+$", ids)) > 0.5

# org.Hs.eg.db bereit?
if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
    !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("R-Pakete für Mapping fehlen.\nInstalliere im Terminal:\n",
       "  conda install -n hipathia-r -c bioconda -c conda-forge ",
       "bioconductor-org.hs.eg.db r-annotationdbi")
}
suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Hs.eg.db) })

# Entrez-Universum und grober Entrez-Check
entrez_universe <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID")
overlap_frac <- mean(ids %in% entrez_universe, na.rm = TRUE)
is_entrez <- looks_numeric && overlap_frac >= 0.5  # konservativ

message(sprintf("ID-Check: looks_numeric=%s, overlap_with_org.Hs.eg.db=%.3f -> is_entrez=%s",
                looks_numeric, overlap_frac, is_entrez))

# 3.2: Robustes Mapping auf Entrez
if (!is_entrez) {
  is_ens    <- grepl("^ENSG\\d+", ids, ignore.case = TRUE)
  is_entrez <- grepl("^[0-9]+$", ids)
  is_symbol <- !(is_ens | is_entrez)

  mapped <- rep(NA_character_, length(ids))

  # SYMBOL -> Entrez
  if (any(is_symbol)) {
    sym_keys <- unique(ids[is_symbol])
    sym2ent  <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys = sym_keys,
                                      keytype = "SYMBOL",
                                      column = "ENTREZID",
                                      multiVals = "first")
    mapped[is_symbol] <- sym2ent[ ids[is_symbol] ]
    message(sprintf("SYMBOL->Entrez: %d/%d zugeordnet",
                    sum(!is.na(mapped[is_symbol])), sum(is_symbol)))
  }

  # ENSEMBL -> Entrez
  if (any(is_ens)) {
    ens_keys <- unique(ids[is_ens])
    ens2ent  <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys = ens_keys,
                                      keytype = "ENSEMBL",
                                      column = "ENTREZID",
                                      multiVals = "first")
    mapped[is_ens] <- ens2ent[ ids[is_ens] ]
    message(sprintf("ENSEMBL->Entrez: %d/%d zugeordnet",
                    sum(!is.na(mapped[is_ens])), sum(is_ens)))
  }

  # bereits Entrez
  mapped[is_entrez] <- ids[is_entrez]

  keep <- !is.na(mapped)
  if (sum(keep) == 0) stop("Kein einziges Gen konnte auf Entrez gemappt werden.")

  expr_mapped <- expr_mat0[keep, , drop = FALSE]
  rownames(expr_mapped) <- mapped[keep]

  # Duplikate (gleiche Entrez) entfernen – erster Eintrag gewinnt
  dup <- duplicated(rownames(expr_mapped))
  if (any(dup)) {
    expr_mapped <- expr_mapped[!dup, , drop = FALSE]
    message(sprintf("Duplikate entfernt: %d", sum(dup)))
  }

  message(sprintf("Nach Mapping: %d Gene (von %d)", nrow(expr_mapped), nrow(expr_mat0)))
} else {
  expr_mapped <- expr_mat0
  message("IDs sind valide Entrez-IDs – kein Mapping nötig.")
}

# 3.3: translate_data (hipathia) + Normalisierung
tr <- hipathia::translate_data(expr_mapped, species = species)
trans_data <- grab_mat(tr)
if (is.null(trans_data) || nrow(trans_data) == 0) {
  stop("translate_data() lieferte 0 Zeilen – passen Organismus/IDs?")
}
message(sprintf("Hipathia-übersetzte Gene: %d", nrow(trans_data)))

exp_data <- hipathia::normalize_data(
  trans_data,
  by_quantiles         = by_quantiles,
  percentil            = percentil_mode,
  truncation_percentil = truncation_percentil
)
message(sprintf("exp_data: %d x %d", nrow(exp_data), ncol(exp_data)))

# 4) LOAD PATHWAYS
message("\n[Load pathways] …")
if (is.null(pathways_list)) {
  pathways <- hipathia::load_pathways(species = species)
} else {
  pathways <- hipathia::load_pathways(species = species, pathways_list = pathways_list)
}
message("Loaded pathways: ", length(hipathia::get_pathways_list(pathways)))

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
results <- hipathia::hipathia(exp_data, pathways, decompose = FALSE)
print(results)

# =====================
# 7) EXTRACT FEATURES (SUBPATHS / FUNCTIONS)
# =====================
path_vals <- hipathia::get_paths_data(results, matrix = TRUE)  # Matrix: subpaths x samples
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

comp <- hipathia::do_wilcoxon(path_vals, sample_group, g1 = group1, g2 = group2)
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

colors_de <- hipathia::node_color_per_de(results, pathways, sample_group, group1, group2)

plot_id <- if (!is.null(pathways_list)) pathways_list[1] else hipathia::get_pathways_list(pathways)[1]
try({
  hipathia::pathway_comparison_plot(comp, metaginfo = pathways, pathway = plot_id, node_colors = colors_de)
}, silent = TRUE)

# =====================
# 12) CREATE REPORT & SAVE SUMMARIES
# =====================
message("\n[Report] Writing HTML report into ", normalizePath(output_dir, mustWork = FALSE))

pathways_summary <- hipathia::get_pathways_summary(comp, pathways)
pathways_summary$percent_significant_paths <- 100 * pathways_summary$num_significant_paths / pmax(1, pathways_summary$num_total_paths)
pathways_summary <- pathways_summary[order(-pathways_summary$percent_significant_paths, -pathways_summary$num_significant_paths), ]
utils::write.table(pathways_summary,
                   file = file.path(output_dir, "pathways_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

report <- hipathia::create_report(comp, pathways, output_dir, node_colors = colors_de)

if (isTRUE(run_local_server)) {
  hipathia::visualize_report(report, port = server_port)
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
