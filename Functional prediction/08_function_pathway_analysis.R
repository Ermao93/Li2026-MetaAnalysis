###############################################################################
# Fish 16S meta-analysis: PICRUSt2 MetaCyc + Tax4Fun2 KEGG Pathway analysis
# Run on cluster with: Rscript 08_function_pathway_analysis.R
#
# Analysis:
#   log10(relative abundance (%) + pseudocount) ~ Host_environment + Region
#   Global Host_environment effect: nested-model ANOVA
#   BH-FDR correction; emmeans adjusted means; pairwise comparisons
#   Top 15 abundance heatmap + up to Top 15 globally significant pathways
#   Pairwise BH-FDR stars are shown inside cells of both heatmaps
###############################################################################

options(stringsAsFactors = FALSE)
set.seed(123)

# ----------------------------- 1. User settings ------------------------------
PROJECT_DIR <- "/public/home/liruojing/fish_16S_merge/analysis_A_global/eukaryota_filtered"
META_FILE <- file.path(PROJECT_DIR, "metadata_diet_region_matched.tsv")
FUNCTION_DIR <- file.path(PROJECT_DIR, "functional_prediction")
OUTPUT_DIR <- file.path(PROJECT_DIR, "08_function_pathway")

GROUP_COL <- "Diet"
REGION_COL <- "Region"
SAMPLE_COL_CANDIDATES <- c("sample-id", "sample_id", "SampleID", "Sample", "sample")

TOP_N_ABUNDANCE <- 15
TOP_N_SIGNIFICANT <- 15
FDR_CUTOFF <- 0.05
PSEUDOCOUNT <- 1e-6

PICRUST_FILE <- file.path(FUNCTION_DIR, "PICRUSt2", "pathways_out",
                          "path_abun_unstrat_descrip.tsv.gz")
TAX4FUN_FILE <- file.path(FUNCTION_DIR, "Tax4Fun2", "run_20261113_215820",
                          "Tax4Fun2_pathway.tsv")
KEGG_INFO_FILE <- "/public/home/liruojing/miniconda3/envs/picrust2_test/lib/python3.10/site-packages/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz"

# ----------------------------- 2. Packages -----------------------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "tibble", "emmeans", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "),
       "。请先在当前 R 环境安装。")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(emmeans)
})

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- 3. Metadata -----------------------------------
if (!file.exists(META_FILE)) stop("找不到 metadata 文件：", META_FILE)
meta <- read.delim(META_FILE, header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
sample_col <- intersect(SAMPLE_COL_CANDIDATES, names(meta))[1]
if (is.na(sample_col)) {
  stop("metadata 中找不到样本 ID 列。当前列名：", paste(names(meta), collapse = " | "))
}
if (!GROUP_COL %in% names(meta)) stop("metadata 缺少分组列：", GROUP_COL)
if (!REGION_COL %in% names(meta)) stop("metadata 缺少区域列：", REGION_COL)

meta[[sample_col]] <- trimws(as.character(meta[[sample_col]]))
meta[[GROUP_COL]] <- droplevels(factor(meta[[GROUP_COL]]))
meta[[REGION_COL]] <- droplevels(factor(meta[[REGION_COL]]))
meta <- meta[!is.na(meta[[sample_col]]) & nzchar(meta[[sample_col]]) &
               !is.na(meta[[GROUP_COL]]) & !is.na(meta[[REGION_COL]]), , drop = FALSE]
meta <- meta[!duplicated(meta[[sample_col]]), , drop = FALSE]

if (nlevels(meta[[GROUP_COL]]) < 2) stop("Host_environment 分组少于 2 组，无法进行组间比较。")
if (nlevels(meta[[REGION_COL]]) < 2) {
  warning("Region 只有一个水平，模型将自动不纳入 Region 项。")
  USE_REGION <- FALSE
} else {
  USE_REGION <- TRUE
}
cat("Metadata 样本数：", nrow(meta), "\n")
cat("Host_environment 分组：", paste(levels(meta[[GROUP_COL]]), collapse = ", "), "\n")
cat("Region 分组：", paste(levels(meta[[REGION_COL]]), collapse = ", "), "\n")

# ----------------------------- 4. Input readers -------------------------------

read_pathway_file <- function(file, meta_samples, id_candidates, description_candidates = character()) {
  if (!file.exists(file)) stop("功能预测文件不存在：", file)
  con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "rt") else base::file(file, "rt")
  on.exit(close(con), add = TRUE)
  dat <- read.delim(con, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (!ncol(dat)) stop("文件没有数据列：", file)

  id_idx <- which(names(dat) %in% id_candidates)[1]
  if (is.na(id_idx)) id_idx <- 1L
  desc_idx <- which(names(dat) %in% description_candidates)[1]
  if (length(desc_idx) == 0 || is.na(desc_idx)) desc_idx <- NA_integer_

  ids <- trimws(as.character(dat[[id_idx]]))
  ids[is.na(ids) | !nzchar(ids)] <- paste0("Feature_", which(is.na(ids) | !nzchar(ids)))
  ids <- make.unique(ids)

  sample_idx <- which(names(dat) %in% meta_samples)
  sample_idx <- setdiff(sample_idx, c(id_idx, desc_idx))
  sample_cols <- names(dat)[sample_idx]
  if (length(sample_cols) < 3) stop("功能表匹配到的 metadata 样本少于 3 个：", file, "\n匹配数：", length(sample_cols))

  mat <- as.matrix(dat[, sample_idx, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  mat[!is.finite(mat) | is.na(mat)] <- 0
  rownames(mat) <- ids

  name_map <- NULL
  if (!is.na(desc_idx)) {
    nm <- trimws(as.character(dat[[desc_idx]]))
    names(nm) <- ids
    name_map <- nm
  }

  list(mat = t(mat), name_map = name_map, id_col = names(dat)[id_idx], description_col = if (!is.na(desc_idx)) names(dat)[desc_idx] else NA_character_, matched_samples = sample_cols)
}


# ----------------------------- 5. KEGG names ----------------------------------
# Use the complete local KEGG pathway description file bundled with PICRUSt2.
# File format: ko00010<space>Glycolysis / Gluconeogenesis
read_local_kegg_names <- function(file) {
  if (!file.exists(file)) stop("找不到本地 KEGG 通路名称文件：", file)
  con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "rt") else base::file(file, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  ids <- sub("[[:space:]].*$", "", lines)
  nms <- sub("^[^[:space:]]+[[:space:]]+", "", lines)
  ok <- grepl("^ko[0-9]{5}$", ids) & nzchar(nms) & nms != lines
  out <- nms[ok]
  names(out) <- ids[ok]
  out <- out[!duplicated(names(out))]
  if (length(out) < 100) stop("本地 KEGG 名称文件解析异常，仅解析到 ", length(out), " 条通路。")
  message("本地 KEGG 通路名称加载成功：", length(out), " 条。")
  out
}

kegg_names <- read_local_kegg_names(KEGG_INFO_FILE)

# ----------------------------- 6. Analysis ------------------------------------
analyse_pathways <- function(input, label, name_map = NULL) {
  mat <- input$mat
  common_samples <- intersect(meta[[sample_col]], rownames(mat))
  if (length(common_samples) < 3) {
    warning(label, " 与 metadata 匹配样本少于 3 个，跳过。")
    return(invisible(NULL))
  }
  md <- meta[match(common_samples, meta[[sample_col]]), , drop = FALSE]
  mm <- mat[common_samples, , drop = FALSE]
  storage.mode(mm) <- "numeric"
  mm[!is.finite(mm) | is.na(mm)] <- 0
  mm[mm < 0] <- 0

  # Drop zero-library samples before relative abundance calculation.
  keep_samples <- rowSums(mm) > 0
  mm <- mm[keep_samples, , drop = FALSE]
  md <- md[keep_samples, , drop = FALSE]
  if (nrow(mm) < 3) {
    warning(label, " 去除零文库样本后不足 3 个，跳过。")
    return(invisible(NULL))
  }

  # Resolve display names; never display pathway IDs on heatmap axes.
  if (is.null(name_map)) {
    pathway_names <- input$name_map[colnames(mm)]
  } else {
    pathway_names <- unname(name_map[colnames(mm)])
  }

  if (is.null(pathway_names)) pathway_names <- rep(NA_character_, ncol(mm))

  pathway_names <- trimws(as.character(pathway_names))
  names(pathway_names) <- colnames(mm)

  # Preserve all pathways in statistical analysis.
  # Exclude unknown or unannotated pathways only from heatmap selection.
  missing_name <- is.na(pathway_names) | !nzchar(pathway_names)
  pathway_names[missing_name] <- paste0("Unannotated_", names(pathway_names)[missing_name])

  valid_heatmap_name <- !grepl(
    "^(unknown|unannotated)(_|$)",
    pathway_names,
    ignore.case = TRUE
  )

  valid_heatmap_ids <- names(pathway_names)[valid_heatmap_name]
  name_coverage <- mean(valid_heatmap_name)

  filter_log <- data.frame(
    Function = names(pathway_names),
    Name = pathway_names,
    Status = ifelse(
      valid_heatmap_name,
      "Retained_annotated",
      "Excluded_from_heatmap_but_retained_in_statistics"
    ),
    Reason = ifelse(
      valid_heatmap_name,
      "",
      "Missing or unknown pathway description; retained in statistical analysis only"
    )
  )

  write.table(
    filter_log,
    file.path(OUTPUT_DIR, paste0(label, "_pathway_annotation_log.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  cat(
    label, " 通路名称覆盖：", sum(valid_heatmap_name), "/", length(valid_heatmap_name),
    "；仅参与统计、不进入热图：", sum(!valid_heatmap_name), "。\n", sep = ""
  )

  pathway_names <- make.unique(pathway_names, sep = " ")
  names(pathway_names) <- colnames(mm)

  # Relative abundance (%).
  rel <- sweep(mm, 1, rowSums(mm), "/") * 100
  rel[!is.finite(rel)] <- 0
  colnames(rel) <- colnames(mm)
  
  md[[GROUP_COL]] <- droplevels(factor(md[[GROUP_COL]]))
  md[[REGION_COL]] <- droplevels(factor(md[[REGION_COL]]))
  if (nlevels(md[[GROUP_COL]]) < 2) stop(label, " 匹配样本中的 Host_environment 少于 2 组。")
  use_region_local <- nlevels(md[[REGION_COL]]) >= 2

  full_terms <- c(GROUP_COL, if (use_region_local) REGION_COL)
  full_formula <- reformulate(full_terms, response = "y")
  reduced_formula <- if (use_region_local) reformulate(REGION_COL, response = "y") else y ~ 1
  design <- data.frame(y = rep(0, nrow(md)))
  design[[GROUP_COL]] <- md[[GROUP_COL]]
  if (use_region_local) design[[REGION_COL]] <- md[[REGION_COL]]
  expected_rank <- qr(model.matrix(full_formula, data = design))$rank

  feature_ids <- colnames(rel)
  global_results <- vector("list", length(feature_ids))
  names(global_results) <- feature_ids
  failed <- character(0)

  cat("\n[", label, "] 开始全局模型检验；样本数=", nrow(rel),
      "；通路数=", ncol(rel), "\n", sep = "")
  for (i in seq_along(feature_ids)) {
    fid <- feature_ids[i]
    y <- log10(rel[, fid] + PSEUDOCOUNT)
    dat <- data.frame(y = y)
    dat[[GROUP_COL]] <- md[[GROUP_COL]]
    if (use_region_local) dat[[REGION_COL]] <- md[[REGION_COL]]
    full <- tryCatch(lm(full_formula, data = dat), error = function(e) NULL)
    reduced <- tryCatch(lm(reduced_formula, data = dat), error = function(e) NULL)
    if (is.null(full) || is.null(reduced) || full$rank < expected_rank ||
        df.residual(full) <= 0) {
      failed[fid] <- "模型拟合失败、秩亏或残差自由度不足"
      next
    }
    aa <- tryCatch(anova(reduced, full), error = function(e) NULL)
    if (is.null(aa) || nrow(aa) < 2 || !is.finite(aa$`Pr(>F)`[2])) {
      failed[fid] <- "嵌套模型比较失败"
      next
    }
    global_results[[fid]] <- data.frame(
      Function = fid, Name = unname(pathway_names[fid]),
      F_value = aa$F[2], P_value = aa$`Pr(>F)`[2]
    )
    if (i %% 500 == 0) cat("  ", i, "/", length(feature_ids), " pathways tested\n", sep = "")
  }

  global_df <- bind_rows(global_results)
  if (!nrow(global_df)) {
    warning(label, " 没有可用的全局检验结果。")
    return(invisible(NULL))
  }
  global_df$FDR <- p.adjust(global_df$P_value, method = "BH")
  global_df <- global_df %>% arrange(FDR, P_value)
  write.table(global_df, file.path(OUTPUT_DIR, paste0(label, "_global_test_BH.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Select the Top N pathways directly by mean relative abundance.
  # PICRUSt2 and Tax4Fun2 are ranked independently; no Top25/manual-review stage.
  mean_rel <- colMeans(rel, na.rm = TRUE)
  abundance_order <- names(sort(mean_rel, decreasing = TRUE))
  abundance_order <- abundance_order[abundance_order %in% valid_heatmap_ids]
  top_ids <- head(abundance_order, TOP_N_ABUNDANCE)
  abundance_top_df <- data.frame(
    Rank = seq_along(top_ids),
    Function = top_ids,
    Name = unname(pathway_names[top_ids]),
    Mean_relative_abundance_percent = as.numeric(mean_rel[top_ids]),
    stringsAsFactors = FALSE
  )
  write.table(abundance_top_df,
              file.path(OUTPUT_DIR, paste0(label, "_Top", length(top_ids), "_abundance_pathways.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # emmeans and pairwise tests are computed for all tested pathways so stars can
  # be drawn on both Top-abundance and significant heatmaps.
  adjusted_list <- list()
  pairwise_list <- list()
  group_levels <- levels(md[[GROUP_COL]])

  cat("[", label, "] 开始 emmeans 和两两比较。\n", sep = "")
  for (i in seq_len(nrow(global_df))) {
    fid <- global_df$Function[i]
    y <- log10(rel[, fid] + PSEUDOCOUNT)
    dat <- data.frame(y = y)
    dat[[GROUP_COL]] <- md[[GROUP_COL]]
    if (use_region_local) dat[[REGION_COL]] <- md[[REGION_COL]]
    fit <- tryCatch(lm(full_formula, data = dat), error = function(e) NULL)
    if (is.null(fit)) next
    emm <- tryCatch(emmeans(fit, specs = GROUP_COL), error = function(e) NULL)
    if (is.null(emm)) next
    ed <- as.data.frame(emm)
    adjusted_list[[fid]] <- data.frame(
      Function = fid, Name = unname(pathway_names[fid]),
      Group = as.character(ed[[GROUP_COL]]),
      Adjusted_mean_log10 = ed$emmean, SE = ed$SE,
      Lower_CL = ed$lower.CL, Upper_CL = ed$upper.CL
    )
    pw <- tryCatch(as.data.frame(pairs(emm, adjust = "none")), error = function(e) NULL)
    if (!is.null(pw) && nrow(pw)) {
      # emmeans contrast labels are "Group1 - Group2"; validate against real levels.
      parts <- strsplit(as.character(pw$contrast), " - ", fixed = TRUE)
      g1 <- vapply(parts, function(z) if (length(z) == 2) z[1] else NA_character_, character(1))
      g2 <- vapply(parts, function(z) if (length(z) == 2) z[2] else NA_character_, character(1))
      valid <- !is.na(g1) & !is.na(g2) & g1 %in% group_levels & g2 %in% group_levels
      if (any(valid)) {
        pairwise_list[[fid]] <- data.frame(
          Function = fid, Name = unname(pathway_names[fid]),
          Group1 = g1[valid], Group2 = g2[valid],
          Contrast = as.character(pw$contrast[valid]),
          Estimate = pw$estimate[valid], SE = pw$SE[valid],
          P_value = pw$p.value[valid]
        )
      }
    }
    if (i %% 200 == 0) cat("  ", i, "/", nrow(global_df), " pathways processed\n", sep = "")
  }

  adjusted_df <- bind_rows(adjusted_list)
  pairwise_df <- bind_rows(pairwise_list)
  if (nrow(adjusted_df)) {
    write.table(adjusted_df, file.path(OUTPUT_DIR, paste0(label, "_adjusted_means.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (nrow(pairwise_df)) {
    # Keep the same global BH adjustment across all pathway x pairwise rows.
    pairwise_df$FDR <- p.adjust(pairwise_df$P_value, method = "BH")
    pairwise_df <- pairwise_df %>% arrange(FDR, P_value)
    write.table(pairwise_df, file.path(OUTPUT_DIR, paste0(label, "_pairwise_BH.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    warning(label, " 没有成功生成两两比较结果；热图不会显示星号。")
  }

  if (length(failed)) {
    failed_df <- data.frame(Function = names(failed), Reason = unname(failed))
    write.table(failed_df, file.path(OUTPUT_DIR, paste0(label, "_model_failed.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Matrix of emmeans adjusted means; rows = pathways, columns = Diet groups.
  make_heatmap_matrix <- function(ids) {
    dd <- adjusted_df %>% filter(Function %in% ids) %>%
      select(Function, Group, Adjusted_mean_log10) %>%
      distinct(Function, Group, .keep_all = TRUE)
    if (!nrow(dd)) return(NULL)
    wide <- tidyr::pivot_wider(dd, names_from = Group, values_from = Adjusted_mean_log10)
    mdf <- as.data.frame(wide)
    rownames(mdf) <- mdf$Function
    mdf$Function <- NULL
    m <- as.matrix(mdf)
    storage.mode(m) <- "numeric"
    m <- m[, group_levels[group_levels %in% colnames(m)], drop = FALSE]
    keep <- apply(m, 1, function(z) all(is.finite(z)) && length(z) >= 2 && sd(z) > 0)
    m <- m[keep, , drop = FALSE]
    if (!nrow(m)) return(NULL)
    z <- t(scale(t(m)))
    z[!is.finite(z)] <- 0
    z
  }

  plot_heatmap <- function(zmat, title, filename, show_pairwise_stars = FALSE) {
    if (is.null(zmat) || !nrow(zmat)) return(invisible(NULL))
    plot_df <- as.data.frame(zmat) %>%
      rownames_to_column("Function") %>%
      pivot_longer(-Function, names_to = "Group", values_to = "Zscore") %>%
      mutate(Name = unname(pathway_names[Function]))
    plot_df$Name <- factor(plot_df$Name, levels = rev(unique(plot_df$Name)))

    if (show_pairwise_stars && nrow(pairwise_df)) {
      stars <- pairwise_df %>%
        filter(Function %in% rownames(zmat), !is.na(FDR), FDR < FDR_CUTOFF) %>%
        mutate(StarLevel = case_when(FDR < 0.001 ~ 3L, FDR < 0.01 ~ 2L, TRUE ~ 1L)) %>%
        select(Function, Group1, Group2, StarLevel) %>%
        pivot_longer(cols = c(Group1, Group2), names_to = "ComparedSide", values_to = "Group") %>%
        group_by(Function, Group) %>%
        summarise(StarLevel = max(StarLevel), .groups = "drop") %>%
        mutate(Star = case_when(StarLevel == 3L ~ "***",
                                StarLevel == 2L ~ "**", TRUE ~ "*"))
      plot_df <- plot_df %>% left_join(stars %>% select(Function, Group, Star),
                                       by = c("Function", "Group"))
    } else {
      plot_df$Star <- ""
    }
    plot_df$Star[is.na(plot_df$Star)] <- ""

    p <- ggplot(plot_df, aes(x = Group, y = Name, fill = Zscore)) +
      geom_tile(colour = "white", linewidth = 0.35) +
      {if (show_pairwise_stars) geom_text(aes(label = Star), colour = "black", size = 3.1, fontface = "bold")} +
      scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
                           midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                           name = "Row Z-score") +
      labs(title = title, x = GROUP_COL, y = NULL) +
      theme_bw(base_size = 10) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(size = 7, colour = "black"),
            axis.title = element_text(colour = "black"),
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8))
    height <- max(6, min(16, 0.30 * nrow(zmat) + 3))
    ggsave(file.path(OUTPUT_DIR, filename), p, width = 9, height = height,
           units = "in", device = "pdf")
    ggsave(file.path(OUTPUT_DIR, sub("\\.pdf$", ".png", filename)), p,
           width = 9, height = height, units = "in", dpi = 300)
    p
  }

  top_z <- make_heatmap_matrix(top_ids)
  plot_heatmap(top_z, paste0(label, " - Top ", TOP_N_ABUNDANCE, " pathways by mean relative abundance"),
               paste0(label, "_Top", length(top_ids), "_abundance_heatmap.pdf"), show_pairwise_stars = TRUE)

sig_ids <- global_df %>%
  filter(!is.na(FDR), FDR < FDR_CUTOFF, Function %in% valid_heatmap_ids) %>%
  arrange(FDR, P_value) %>%
  slice_head(n = TOP_N_SIGNIFICANT) %>%
  pull(Function)
  
  if (length(sig_ids)) {
    sig_z <- make_heatmap_matrix(sig_ids)
    plot_heatmap(sig_z, paste0(label, " - Top ", length(sig_ids), " significant pathways (global BH-FDR < ", FDR_CUTOFF, ")"),
                 paste0(label, "_Top", length(sig_ids), "_significant_heatmap.pdf"), show_pairwise_stars = TRUE)
  } else {
    message(label, "：没有全局 BH-FDR < ", FDR_CUTOFF, " 的通路，不生成显著差异热图。")
  }

  summary_df <- data.frame(
    Analysis = label, Matched_samples = nrow(rel),
    Matched_samples = ncol(rel), Tested_pathways = nrow(global_df),
    Global_significant = sum(global_df$FDR < FDR_CUTOFF, na.rm = TRUE),
    Pairwise_comparisons = nrow(pairwise_df),
    KEGG_or_MetaCyc_name_coverage = round(name_coverage * 100, 2)
  )
  write.table(summary_df, file.path(OUTPUT_DIR, paste0(label, "_summary.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("[", label, "] 完成。结果目录：", OUTPUT_DIR, "\n", sep = "")
  invisible(list(global = global_df, adjusted = adjusted_df, pairwise = pairwise_df))
}

# ----------------------------- 7. Run both ------------------------------------
cat("\n========== PICRUSt2 MetaCyc Pathway ==========\n")
picrust <- read_pathway_file(
  PICRUST_FILE, meta[[sample_col]],
  id_candidates = c("pathway", "Pathway", "feature", "Feature_ID"),
  description_candidates = c("description", "Description", "name", "Name")
)
if (is.null(picrust$name_map)) stop("PICRUSt2 描述文件中未检测到 description 列。")
cat("PICRUSt2 匹配样本数：", length(picrust$matched_samples), "\n")
cat("PICRUSt2 通路描述覆盖：",
    sum(!is.na(picrust$name_map) & nzchar(picrust$name_map)), "/", length(picrust$name_map), "\n")
picrust_result <- analyse_pathways(picrust, "PICRUSt2_MetaCyc_Pathway")

cat("\n========== Tax4Fun2 KEGG Pathway ==========\n")
tax4fun <- read_pathway_file(
  TAX4FUN_FILE, meta[[sample_col]],
  id_candidates = c("pathway", "Pathway", "ko", "KO", "Feature_ID"),
  description_candidates = c("description", "Description", "name", "Name")
)
# Map each observed Tax4Fun2 pathway ID to the complete local KEGG dictionary.
tax_ids <- colnames(tax4fun$mat)
tax_names <- unname(kegg_names[tax_ids])
names(tax_names) <- tax_ids
mapped <- !is.na(tax_names) & nzchar(tax_names)
cat("Tax4Fun2 匹配样本数：", length(tax4fun$matched_samples), "\n")
cat("Tax4Fun2 KEGG 名称覆盖：", sum(mapped), "/", length(tax_ids),
    " (", round(100 * mean(mapped), 2), "%)\n", sep = "")
# Export mapping audit; unmapped IDs remain in the analysis and are labelled by ID.
tax_map_audit <- data.frame(
  Function = tax_ids,
  Name = ifelse(mapped, tax_names, paste0("Unannotated_", tax_ids)),
  Mapped = mapped,
  stringsAsFactors = FALSE
)
write.table(tax_map_audit,
            file.path(OUTPUT_DIR, "Tax4Fun2_KEGG_pathway_name_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (any(!mapped)) {
  write.table(tax_map_audit[!mapped, , drop = FALSE],
              file.path(OUTPUT_DIR, "Tax4Fun2_unmapped_pathways.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
tax4fun$name_map <- tax_names
tax4fun_result <- analyse_pathways(tax4fun, "Tax4Fun2_KEGG_Pathway",
                                   name_map = tax4fun$name_map)

cat("\n全部功能预测 Pathway 分析完成。\n结果目录：", OUTPUT_DIR, "\n", sep = "")
