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
pathways_list <- NULL
            # e.g., c("hsa03320","hsa04014") or NULL for all

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
    stop("Fehlendes Paket: ", pkg, ". Bitte im Environment installieren (Conda/BiocManager).")
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
# --- ID-Diagnose ---
sym_like   <- mean(grepl("^[A-Za-z]", rownames(expr_mat0))) > 0.3
ens_like   <- mean(grepl("^ENSG", rownames(expr_mat0))) > 0.3
entrez_like<- mean(grepl("^[0-9]+$", rownames(expr_mat0))) > 0.3
message(sprintf("ID-Heuristik ~ SYMBOL:%s ENSEMBL:%s ENTREZ:%s", sym_like, ens_like, entrez_like))

tic <- function() assign(".tic", Sys.time(), .GlobalEnv)
toc <- function(lbl="step") {
  t <- get(".tic", envir=.GlobalEnv); d <- difftime(Sys.time(), t, units="mins")
  message(sprintf("[%s] %.2f min", lbl, as.numeric(d)))
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

# --- Gruppen-Check ---
stopifnot("group" %in% colnames(col_data))
grp_levels <- unique(as.character(col_data$group))
if (!all(c(group1, group2) %in% grp_levels)) {
  stop(sprintf("Gruppen passen nicht: group1/group2 = %s/%s, vorhandene: %s",
               group1, group2, paste(grp_levels, collapse=", ")))
}

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
# --- Aufräumen: NAs raus/imputieren, konstante Subpaths erkennen ---
# 2.1 Zeilen, die komplett NA sind, entfernen
keep_any <- rowSums(is.finite(path_vals)) > 1
if (sum(keep_any) == 0) {
  stop("Alle Subpaths sind komplett NA. Lade mehr Pfade (pathways_list <- NULL) oder prüfe Mapping/normalize_data.")
}
path_vals <- path_vals[keep_any, , drop = FALSE]

# 2.2 NA zeilenweise mit Zeilenmittelwert imputieren (falls ganze Zeile NA war, wäre sie oben schon raus)
row_means <- rowMeans(path_vals, na.rm = TRUE)
# ersetze NA-Zellen durch den jeweiligen Zeilenmittelwert
na_idx <- which(is.na(path_vals), arr.ind = TRUE)
if (length(na_idx) > 0) {
  path_vals[na_idx] <- row_means[na_idx[,1]]
}

# 2.3 Konstantheit prüfen
vars <- apply(path_vals, 1, var)
nzv  <- vars > 0 & is.finite(vars)
if (sum(nzv) < 2) {
  warning("Zu wenig variable Subpaths für PCA nach NA-Handling. Ich überspringe PCA.")
} else {
  path_vals <- path_vals[nzv, , drop = FALSE]
}
# --- Diagnose-Log nach Cleanup ---
message(sprintf("Subpaths nach Cleanup: %d (von vorher %d)", nrow(path_vals), length(vars)))
message(sprintf("NA-Zellen nach Imputation: %d", sum(is.na(path_vals))))

# --- Option: schwach variable Subpaths filtern (QoL) ---
# Deaktiviert per default; aktiviere bei Bedarf.
filter_low_var <- FALSE
if (isTRUE(filter_low_var)) {
  v <- apply(path_vals, 1, stats::var)
  keep <- v > stats::quantile(v, 0.2, na.rm = TRUE)  # oberste 80% Varianz behalten
  path_vals <- path_vals[keep, , drop = FALSE]
  message(sprintf("Varianzfilter aktiv: %d/%d Subpaths behalten", sum(keep), length(v)))
}

# =====================
# 8) STATS: WILCOXON g1 VS g2
# =====================
message("
[Wilcoxon] ", group1, " vs ", group2)
sample_group <- col_data$group
names(sample_group) <- rownames(col_data)
stopifnot(all(names(sample_group) == colnames(path_vals)))

tic <- function() assign(".tic", Sys.time(), .GlobalEnv)
toc <- function(lbl="step") {
  t <- get(".tic", envir=.GlobalEnv); d <- difftime(Sys.time(), t, units="mins")
  message(sprintf("[%s] %.2f min", lbl, as.numeric(d)))
}

tic(); tr <- hipathia::translate_data(expr_mapped, species = species); toc("translate_data")
tic(); results <- hipathia::hipathia(exp_data, pathways, decompose = FALSE); toc("hipathia")
tic(); comp <- hipathia::do_wilcoxon(path_vals, sample_group, g1 = group1, g2 = group2); toc("wilcoxon")


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

# --- Komfort: Top-Tabellen exportieren ---
ord_p  <- order(comp$p.value, decreasing = FALSE)
ord_q  <- order(comp$q.value, decreasing = FALSE)
utils::write.table(comp[head(ord_p, 50), ],
  file = file.path(output_dir, "subpath_top50_by_p.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE
)
utils::write.table(comp[head(ord_q, 50), ],
  file = file.path(output_dir, "subpath_top50_by_q.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE
)
message("Top-5 p: ", paste(signif(comp$p.value[head(ord_p,5)], 3), collapse=", "))
message("Top-5 q: ", paste(signif(comp$q.value[head(ord_q,5)], 3), collapse=", "))

# =====================
# PCA (robust via prcomp)
# =====================
if (use_pca) {
  message("\n[PCA] (robust via prcomp) …")

  if (nrow(path_vals) < 2 || ncol(path_vals) < 2) {
    warning("Zu wenige Dimensionen für PCA (", nrow(path_vals), "x", ncol(path_vals), "). Überspringe PCA.")
  } else {
    # optional: auf Top-Varianz featuren beschränken
    vars <- apply(path_vals, 1, var)
    k_total <- length(vars)
    k_user  <- if (is.na(max_pca_features)) k_total else min(k_total, max_pca_features)
    top_ids <- names(sort(vars, decreasing = TRUE))[seq_len(k_user)]
    X <- path_vals[top_ids, , drop = FALSE]               # Subpaths x Samples

    # Gruppen-Vektor passend zu den Spalten
    sample_group <- col_data$group
    names(sample_group) <- rownames(col_data)
    grp <- factor(sample_group[colnames(X)], levels = c(group1, group2))

    # PCA auf z-skalierten Daten; t(X) = Samples x Features
    pca_res <- prcomp(t(X), center = TRUE, scale. = TRUE)
    imp <- summary(pca_res)$importance[2, 1:2] * 100

    pdf(file.path(output_dir, "pca_prcomp.pdf"))
    cols <- as.integer(grp)
    plot(pca_res$x[,1], pca_res$x[,2],
         col = cols, pch = 19,
         xlab = sprintf("PC1 (%.1f%%)", imp[1]),
         ylab = sprintf("PC2 (%.1f%%)", imp[2]),
         main = sprintf("PCA on Hipathia subpaths (n=%d, k=%d)", ncol(X), nrow(X)))
    legend("topright", legend = levels(grp), col = seq_along(levels(grp)), pch = 19)
    dev.off()
    message("PCA gespeichert: ", file.path(output_dir, "pca_prcomp.pdf"))
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

# ---- Kompatibilität für create_report() (alte Hipathia-Erwartungen) ----
# 1) FDR-Spalte so benennen, wie der Report sie erwartet
if (!"FDRp.value" %in% colnames(comp)) {
  comp$FDRp.value <- if ("q.value" %in% colnames(comp)) comp$q.value else p.adjust(comp$p.value, method = p_adj_method)
}

# 2) Richtung bestimmen (g1 - g2) und "status" setzen
g1_means <- rowMeans(path_vals[, sample_group == group1, drop = FALSE], na.rm = TRUE)
g2_means <- rowMeans(path_vals[, sample_group == group2, drop = FALSE], na.rm = TRUE)

# auf comp-Rownames ausrichten
comp$mean_g1 <- g1_means[rownames(comp)]
comp$mean_g2 <- g2_means[rownames(comp)]
comp$effect  <- comp$mean_g1 - comp$mean_g2
comp$status  <- ifelse(comp$effect > 0, "UP", "DOWN")

# 3) (optional) NAs in FDR vermeiden
comp$FDRp.value[!is.finite(comp$FDRp.value)] <- 1


report <- hipathia::create_report(comp, pathways, output_dir,
                                  node_colors = colors_de, conf = alpha)

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

# --- Session Info & Repro-Log ---
writeLines(c(
  sprintf("Date: %s", Sys.time()),
  sprintf("Working dir: %s", getwd()),
  sprintf("R version: %s", R.version.string)
), con = file.path(output_dir, "run_info.txt"))
sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
saveRDS(list(
  results = results,
  comp    = comp,
  path_vals = path_vals,
  pathways = pathways,
  col_data = col_data
), file = file.path(output_dir, "objects_for_repro.rds"))

message("\n[DONE]")
