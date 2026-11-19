###############################################################################
# Fish 16S meta-analysis - Part 1 (unified, modular, directly runnable)
# ----------------------------------------------------------------------------
# 使用方法：
#   1) 只需修改下方 "0. USER SETTINGS" 中的 package_dir 与 group_col。
#   2) 按需开关各 run_* 模块。
#   3) 全部结果输出到： <package_dir>/results_part1/<group_col>/  (可用环境变量覆盖)
#
# 与旧版相比的修正 / 改进：
#   * LEfSe 改为自包含实现（Kruskal-Wallis + 组间 CPM 均值 LDA 效应量），不再依赖
#     microbiomeMarker / phyloseq，彻底避免其空 marker_table 的 validObject 报错。
#   * 按真实 PCoA 文件格式(Eigvals / Proportion explained / Species / Site / Biplot
#     / Site constraints)重写 read_pcoa，不再假设旧的表头结构。
#   * 两个 BIOM(稀释表与 otu_table_no_eukaryota)实际均为 HDF5 格式，统一用
#     rhdf5 直接重建稀疏矩阵读取，不再依赖 biomformat 对 HDF5 的不稳定支持。
#   * 稀释表用稀疏矩阵聚合 ASV -> Phylum / Genus，避免把 35686 x 2545 稠密化。
#   * ANCOM-BC2 结果列名按 2.12.x 版本修正(全局: diff_abn；成对: diff_robust_)，
#     且分组系数列名按 group_col 动态识别，换分组变量无需再改。
#   * 功能预测(FAPROTAX / Tax4Fun2 KO / pathway)统一改为 代码3 的协变量校正方法：
#     lm(log(相对丰度) ~ group + Region) 嵌套比较 + emmeans 校正均值 + 两两比较 +
#     BH-FDR + 热图。KO 数量大，仅对全局显著项做两两比较。
#   * 全部图形统一 theme_paper() 与统一 group_palette。
###############################################################################

############################ 0. USER SETTINGS ################################

package_dir <- "F:/湛江湾实验室/文章/meta综述/16s下游数据分析11.3/fish_16S_analysis_package_20261116"

# 核心分组变量。可改为 Diet、WL、BS、MP、Host_environment、LR、SP 等分类变量。
group_col <- "Host_environment"

# 技术扩增区域协变量（PERMANOVA / 功能分析 / ANCOM-BC2 中作为协变量校正）。
region_col <- "Region"

# ANCOM-BC2 / LEfSe 的参考水平（仅当该值确实存在于分组水平中时才生效）。
# 若分组不是 Diet 或不想指定，可设为 NULL 使用默认水平。
reference_level <- "Carnivore"

# 可选：手动指定分组水平顺序；NULL 表示自动(按字母序，参考水平会被置前)。
group_levels_manual <- NULL

# 是否运行各模块
run_alpha                <- TRUE
run_beta                 <- TRUE
run_composition          <- TRUE
run_genus_differential   <- TRUE
run_core                 <- TRUE
run_venn                 <- TRUE
run_lefse                <- TRUE
run_functional_analysis  <- TRUE
run_continuous_traits    <- TRUE
run_ancombc2             <- TRUE

# 参数
permutations     <- 999
top_genus_n      <- 20
top_phylum_n     <- 12
top_diff_n       <- 30
top_function_n   <- 20
fdr_cutoff       <- 0.05
core_prevalence  <- 0.70
core_top_n       <- 50
lefse_kw_cutoff  <- 0.05 # LEfSe 的 Kruskal-Wallis 阈值（决定哪些属计算 LDA 效应量）
lefse_final_lda  <- 4    # 最终差异属后筛选的 |LDA| 阈值（与 KW-BH 双重过滤）
pseudocount      <- 1e-6
recalculate_nmds <- FALSE   # 是否重算 NMDS（默认 FALSE，优先用数据包内已有坐标）

set.seed(123)
options(stringsAsFactors = FALSE)

############################ 1. PACKAGES ######################################

# 按启用模块收集所需包，只安装真正缺失的包；安装失败不中断（离线时警告并继续）。
cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "tibble", "stringr", "vegan",
               "RColorBrewer", "scales")
if (run_venn) cran_pkgs <- c(cran_pkgs, "ggvenn")
if (run_core) cran_pkgs <- c(cran_pkgs, "viridis")
if (run_functional_analysis) cran_pkgs <- c(cran_pkgs, "emmeans")

bioc_pkgs <- c("rhdf5")                            # 两个 BIOM 均为 HDF5 格式，读取始终需要 rhdf5
if (run_ancombc2) bioc_pkgs <- c(bioc_pkgs, "ANCOMBC")

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(TRUE))
  tryCatch({
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(missing, ask = FALSE, update = FALSE)
    } else {
      install.packages(missing, repos = "https://cloud.r-project.org")
    }
  }, error = function(e) {
    warning("以下包安装失败（可能无网络），相关模块将跳过：",
            paste(missing, collapse = ", "), " - ", conditionMessage(e))
  })
  invisible(TRUE)
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(vegan); library(RColorBrewer); library(scales)
  library(Matrix)
})

############################ 2. PATHS AND UTILITIES ############################

if (!dir.exists(package_dir)) stop("找不到数据包目录：", package_dir)

# 输出目录默认放在数据包内；可用环境变量 FISH16S_OUT_DIR 覆盖(便于测试/重定向)。
out_dir <- Sys.getenv("FISH16S_OUT_DIR",
                      unset = file.path(package_dir, "results_part1", group_col))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p   <- function(...) file.path(package_dir, ...)
out <- function(...) file.path(out_dir, ...)
make_dir <- function(...) { z <- out(...); dir.create(z, recursive = TRUE, showWarnings = FALSE); z }

dirs <- c("01_alpha", "02_beta/Figures", "02_beta/Tables", "03_composition",
          "04_genus_differential", "05_core", "06_venn", "07_LEfSe",
          "08_function", "09_continuous_traits", "10_ANCOMBC2",
          "11_tables", "12_QC")
invisible(lapply(dirs, make_dir))

read_tsv <- function(file) {
  read.delim(file, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_table <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

save_plot <- function(plot, filename, width = 8, height = 6) {
  pdf_file <- if (grepl("\\.pdf$", filename, ignore.case = TRUE)) filename else paste0(filename, ".pdf")
  png_file <- sub("\\.pdf$", ".png", pdf_file, ignore.case = TRUE)
  dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(pdf_file, plot, width = width, height = height, device = "pdf", limitsize = FALSE)
  ggsave(png_file, plot, width = width, height = height, dpi = 300, limitsize = FALSE)
}

theme_paper <- function() {
  theme_bw(base_size = 12) +
    theme(text = element_text(colour = "black"),
          axis.text = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          plot.title = element_text(colour = "black", face = "bold", hjust = 0.5),
          plot.subtitle = element_text(colour = "black", hjust = 0.5),
          legend.title = element_text(colour = "black"),
          legend.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"),
          strip.text = element_text(colour = "black", face = "bold"))
}

make_palette <- function(x) {
  x <- unique(as.character(x))
  n <- length(x)
  cols <- if (n <= 8) brewer.pal(max(3, n), "Set1")[seq_len(n)] else
    colorRampPalette(brewer.pal(8, "Set1"))(n)
  names(cols) <- x
  cols
}

kw_p <- function(y, g) {
  ok <- is.finite(y) & !is.na(g)
  y <- y[ok]; g <- droplevels(as.factor(g[ok]))
  if (length(unique(g)) < 2) return(NA_real_)
  tryCatch(stats::kruskal.test(y, g)$p.value, error = function(e) NA_real_)
}

relative_abundance <- function(mat) {
  mat <- as.matrix(mat)
  den <- rowSums(mat); den[den == 0] <- 1
  sweep(mat, 1, den, "/")
}

# 读取 HDF5 格式 BIOM，返回 features x samples 的 dgCMatrix。
# 本数据包的两个 .biom 文件（稀释表与 no_eukaryota 原始表）均为 QIIME2 导出的
# BIOM 2.1 HDF5 格式，统一用此函数直接重建稀疏矩阵。
read_hdf5_biom_sparse <- function(file) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) stop("读取 HDF5 BIOM 需要 rhdf5 包。")
  obs_ids    <- as.character(rhdf5::h5read(file, "/observation/ids"))
  sample_ids <- as.character(rhdf5::h5read(file, "/sample/ids"))
  data       <- as.numeric(rhdf5::h5read(file, "/observation/matrix/data"))
  indices    <- as.integer(rhdf5::h5read(file, "/observation/matrix/indices"))
  indptr     <- as.integer(rhdf5::h5read(file, "/observation/matrix/indptr"))
  n_row <- length(obs_ids); n_col <- length(sample_ids)
  i <- rep.int(seq_len(n_row), diff(indptr))   # CSR: 行=feature
  j <- indices + 1L                            # 列=sample
  M <- Matrix::sparseMatrix(i = i, j = j, x = data, dims = c(n_row, n_col),
                            dimnames = list(obs_ids, sample_ids))
  M
}

############################ 3. INPUT CHECKS AND METADATA #####################

required_files <- c(
  "01_Metadata/metadata_full.tsv",
  "01_Metadata/metadata_rarefied_1000.tsv",
  "02_Abundance/BIOM/otu_table_rarefied_1000.biom",
  "03_Taxonomy_Sequences/taxonomy_no_eukaryota.tsv"
)
if (run_alpha) {
  required_files <- c(required_files,
    paste0("05_Diversity/alpha/", c("shannon", "observed_features", "pielou_e",
                                    "faith_pd", "simpson"), "/alpha-diversity.tsv"))
}
if (run_ancombc2) required_files <- c(required_files, "02_Abundance/BIOM/otu_table_no_eukaryota.biom")
missing <- required_files[!file.exists(p(required_files))]
if (length(missing)) stop("数据包缺少当前启用模块所需文件：\n", paste(missing, collapse = "\n"))

meta_full <- read_tsv(p("01_Metadata/metadata_full.tsv")) %>% distinct(`sample-id`, .keep_all = TRUE)
meta_rare <- read_tsv(p("01_Metadata/metadata_rarefied_1000.tsv")) %>% distinct(`sample-id`, .keep_all = TRUE)

if (!group_col %in% names(meta_rare)) stop("分组变量不存在于 metadata_rarefied_1000.tsv：", group_col)
if ((run_functional_analysis || run_ancombc2) && !group_col %in% names(meta_full))
  stop("分组变量不存在于 metadata_full.tsv：", group_col)

if (is.numeric(meta_rare[[group_col]]))
  stop("group_col 是数值型变量；若它是分类组，请先在 metadata 中转为字符/因子。")

if (!is.null(group_levels_manual)) {
  meta_rare[[group_col]] <- factor(meta_rare[[group_col]], levels = group_levels_manual)
} else {
  meta_rare[[group_col]] <- factor(meta_rare[[group_col]])
}
if (group_col %in% names(meta_full)) meta_full[[group_col]] <- factor(meta_full[[group_col]])
if (!is.null(reference_level) && reference_level %in% levels(meta_rare[[group_col]]))
  meta_rare[[group_col]] <- relevel(meta_rare[[group_col]], ref = reference_level)
if (region_col %in% names(meta_rare)) meta_rare[[region_col]] <- factor(meta_rare[[region_col]])
if (region_col %in% names(meta_full))  meta_full[[region_col]]  <- factor(meta_full[[region_col]])
if (!region_col %in% names(meta_rare)) warning("元数据缺少 Region；相关模型将不纳入 Region。")

meta_rare[[group_col]] <- droplevels(meta_rare[[group_col]])
if (nlevels(meta_rare[[group_col]]) < 2) stop("分组变量至少需要两个有效水平：", group_col)

group_levels <- levels(meta_rare[[group_col]])
group_palette <- make_palette(group_levels)

for (v in intersect(c("TL", "Max_length"), names(meta_rare))) meta_rare[[v]] <- suppressWarnings(as.numeric(meta_rare[[v]]))
for (v in intersect(c("TL", "Max_length"), names(meta_full))) meta_full[[v]] <- suppressWarnings(as.numeric(meta_full[[v]]))

write_table(data.frame(
  item = c("metadata_full_samples", "metadata_rarefied_samples", "group_variable",
           "region_variable", "reference_level", "group_levels"),
  value = c(nrow(meta_full), nrow(meta_rare), group_col, region_col,
            ifelse(is.null(reference_level), "automatic", reference_level),
            paste(group_levels, collapse = ";"))
), out("12_QC/metadata_summary.tsv"))

############################ 4. BIOM + TAXONOMY + GENUS MATRICES ##############

otu0 <- read_hdf5_biom_sparse(p("02_Abundance/BIOM/otu_table_rarefied_1000.biom"))  # feature x sample
if (!all(colnames(otu0) %in% meta_rare[["sample-id"]]))
  stop("稀释 BIOM 样本 ID 无法与 metadata_rarefied_1000.tsv 匹配。")
otu <- Matrix::t(otu0)                              # sample x feature
common_samples <- intersect(meta_rare[["sample-id"]], rownames(otu))
if (length(common_samples) < 3) stop("稀释 BIOM 与元数据匹配样本不足。")
meta_rare <- meta_rare[match(common_samples, meta_rare[["sample-id"]]), , drop = FALSE]
meta_rare[[group_col]] <- droplevels(meta_rare[[group_col]])
group_levels <- levels(meta_rare[[group_col]])
group_palette <- make_palette(group_levels)
otu <- otu[common_samples, , drop = FALSE]
otu <- otu[, Matrix::colSums(otu) > 0, drop = FALSE]

tax <- read_tsv(p("03_Taxonomy_Sequences/taxonomy_no_eukaryota.tsv"))
feature_col <- intersect(c("Feature ID", "Feature.ID", "feature-id", "FeatureID", "ASV"), names(tax))[1]
if (is.na(feature_col)) feature_col <- names(tax)[1]
taxon_col <- intersect(c("Taxon", "taxonomy", "Taxonomy"), names(tax))[1]
if (is.na(taxon_col)) stop("taxonomy 表中找不到 Taxon 分类字符串列。")

tax2 <- tax %>%
  transmute(FeatureID = as.character(.data[[feature_col]]),
            TaxonString = as.character(.data[[taxon_col]])) %>%
  distinct(FeatureID, .keep_all = TRUE)

split_tax <- function(x) if (is.na(x) || !nzchar(x)) character(0) else strsplit(x, "\\s*[;|]\\s*")[[1]]
extract_rank <- function(x, prefix, index) {
  z <- split_tax(x)
  hit <- z[grepl(paste0("^", prefix), z, ignore.case = TRUE)]
  ans <- if (length(hit)) sub(paste0("^", prefix), "", tail(hit, 1), ignore.case = TRUE)
         else if (length(z) >= index) z[index] else NA_character_
  ans <- trimws(ans)
  if (is.na(ans) || !nzchar(ans) ||
      grepl("unclassified|uncultured|unknown|metagenome|unassigned|incertae_sedis", ans, ignore.case = TRUE))
    "Unclassified" else ans
}
tax2$Phylum <- vapply(tax2$TaxonString, extract_rank, character(1), prefix = "p__", index = 2)
tax2$Genus  <- vapply(tax2$TaxonString, extract_rank, character(1), prefix = "g__", index = 6)

feature_tax <- tax2 %>% filter(FeatureID %in% colnames(otu))
if (!nrow(feature_tax)) stop("taxonomy 与稀释 BIOM 没有匹配 Feature ID。")
keep_feat <- colnames(otu) %in% feature_tax$FeatureID
otu_tax <- otu[, keep_feat, drop = FALSE]

aggregate_sparse <- function(mat, labels) {
  labels[is.na(labels) | !nzchar(labels)] <- "Unclassified"
  u <- sort(unique(labels))
  map <- Matrix::sparseMatrix(i = match(labels, u), j = seq_along(labels), x = 1,
                              dims = c(length(u), length(labels)))
  res <- mat %*% Matrix::t(map)
  colnames(res) <- u
  as.matrix(res)
}

phylum_lab <- feature_tax$Phylum[match(colnames(otu_tax), feature_tax$FeatureID)]
genus_lab  <- feature_tax$Genus[match(colnames(otu_tax), feature_tax$FeatureID)]
phylum_mat <- aggregate_sparse(otu_tax, phylum_lab)
genus_mat  <- aggregate_sparse(otu_tax, genus_lab)
genus_mat  <- genus_mat[, colSums(genus_mat) > 0, drop = FALSE]

write_table(data.frame(
  item = c("matched_samples", "retained_ASVs", "taxonomy_matched_ASVs",
           "phylum_number", "genus_number"),
  value = c(nrow(otu), ncol(otu), ncol(otu_tax), ncol(phylum_mat), ncol(genus_mat))
), out("12_QC/rarefied_BIOM_and_taxonomy_summary.tsv"))

############################ 5. ALPHA DIVERSITY ################################

alpha_long <- data.frame(metric = character(), value = numeric())
if (run_alpha) {
  alpha_names <- c(shannon = "Shannon", observed_features = "Observed_features",
                   pielou_e = "Pielou_evenness", faith_pd = "Faith_PD",
                   simpson = "Simpson")
  read_alpha <- function(metric, label) {
    x <- read_tsv(p("05_Diversity", "alpha", metric, "alpha-diversity.tsv"))
    names(x)[1] <- "sample-id"
    data.frame(`sample-id` = as.character(x[[1]]),
               value = as.numeric(x[[2]]),
               metric = label, check.names = FALSE)
  }
  alpha_long <- bind_rows(Map(read_alpha, names(alpha_names), alpha_names)) %>%
    left_join(meta_rare, by = "sample-id") %>%
    filter(!is.na(.data[[group_col]]), is.finite(value))

  alpha_stats <- alpha_long %>%
    group_by(metric) %>%
    summarise(n = n(), p_value = kw_p(value, .data[[group_col]]), .groups = "drop") %>%
    mutate(FDR = p.adjust(p_value, "BH"))

  write_table(alpha_stats, out("11_tables/alpha_KW_BH.tsv"))
  write_table(alpha_long, out("11_tables/alpha_diversity_all_samples.tsv"))

  for (m in unique(alpha_long$metric)) {
    dat <- filter(alpha_long, metric == m)
    st  <- filter(alpha_stats, metric == m)
    g <- ggplot(dat, aes(x = .data[[group_col]], y = value, fill = .data[[group_col]])) +
      geom_boxplot(width = 0.68, outlier.shape = NA, alpha = 0.8) +
      geom_jitter(width = 0.15, size = 0.6, alpha = 0.35) +
      scale_fill_manual(values = group_palette, drop = FALSE) +
      labs(x = group_col, y = m, title = m,
           subtitle = paste0("Kruskal-Wallis P = ", signif(st$p_value, 3),
                             "; BH-FDR = ", signif(st$FDR, 3))) +
      theme_paper() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
    save_plot(g, out("01_alpha", paste0(m, ".pdf")), 7, 5.5)
  }
}

############################ 6. BETA DIVERSITY #################################

distance_names <- c(Bray_Curtis = "bray_curtis", Jaccard = "jaccard",
                    Weighted_UniFrac = "weighted_unifrac",
                    Unweighted_UniFrac = "unweighted_unifrac")

read_distance <- function(metric) {
  f <- p("05_Diversity", "beta", metric, "distance-matrix.tsv")
  if (!file.exists(f)) stop("距离矩阵不存在：", f)
  x <- read_tsv(f)
  ids <- as.character(x[[1]])
  mat <- as.matrix(x[, -1, drop = FALSE]); storage.mode(mat) <- "numeric"
  rownames(mat) <- ids; colnames(mat) <- names(x)[-1]
  common <- intersect(rownames(mat), colnames(mat))
  mat <- mat[common, common, drop = FALSE]
  if (length(common) < 3 || any(!is.finite(mat)) || !identical(rownames(mat), colnames(mat)))
    stop("距离矩阵格式或样本异常：", metric)
  mat
}

read_pcoa <- function(metric) {
  f <- p("05_Diversity", "pcoa", metric, "ordination.txt")
  if (!file.exists(f)) stop("PCoA 文件不存在：", f)
  lines <- readLines(f, warn = FALSE)
  pp <- which(startsWith(lines, "Proportion explained"))
  ss <- which(startsWith(lines, "Site\t"))
  if (length(pp) != 1 || length(ss) != 1) stop("无法识别 PCoA 文件结构：", f)
  prop <- suppressWarnings(as.numeric(strsplit(lines[pp + 1], "\t", fixed = TRUE)[[1]]))
  prop <- prop[is.finite(prop)]
  site_info <- strsplit(lines[ss], "\t", fixed = TRUE)[[1]]
  n_sites <- as.integer(site_info[2])
  if (is.na(n_sites) || n_sites < 3 || length(prop) < 2) stop("PCoA 坐标或解释率不足：", f)
  site_rows <- strsplit(lines[(ss + 1):(ss + n_sites)], "\t", fixed = TRUE)
  site_rows <- Filter(function(z) length(z) >= 3 &&
                        all(!is.na(suppressWarnings(as.numeric(z[2:3])))), site_rows)
  dat <- data.frame(`sample-id` = vapply(site_rows, `[`, character(1), 1),
                    PCoA1 = as.numeric(vapply(site_rows, `[`, character(1), 2)),
                    PCoA2 = as.numeric(vapply(site_rows, `[`, character(1), 3)),
                    check.names = FALSE)
  list(data = dat, prop = prop)
}

read_nmds <- function(metric) {
  f <- p("05_Diversity", "nmds", paste0(metric, "_coordinates.tsv"))
  if (!file.exists(f)) return(NULL)
  x <- read_tsv(f)
  if (!all(c("sample_id", "NMDS1", "NMDS2") %in% names(x))) return(NULL)
  x %>% transmute(`sample-id` = as.character(sample_id),
                  NMDS1 = as.numeric(NMDS1), NMDS2 = as.numeric(NMDS2))
}

if (run_beta) {
  meta_beta <- meta_rare %>% filter(!is.na(.data[[group_col]]))
  if (region_col %in% names(meta_beta)) meta_beta <- meta_beta %>% filter(!is.na(.data[[region_col]]))

  permanova_results <- list(); permdisp_results <- list()
  for (lab in names(distance_names)) {
    metric <- distance_names[[lab]]
    dmat <- tryCatch(read_distance(metric), error = function(e) { warning(conditionMessage(e)); NULL })
    if (is.null(dmat)) next
    ids <- intersect(rownames(dmat), meta_beta[["sample-id"]])
    if (length(ids) < 3) next
    d  <- as.dist(dmat[ids, ids, drop = FALSE])
    md <- meta_beta[match(ids, meta_beta[["sample-id"]]), , drop = FALSE]

    # PCoA
    pc <- tryCatch(read_pcoa(metric), error = function(e) { warning(conditionMessage(e)); NULL })
    if (!is.null(pc)) {
      pd <- inner_join(pc$data, md, by = "sample-id")
      if (nrow(pd) >= 3) {
        g <- ggplot(pd, aes(PCoA1, PCoA2, colour = .data[[group_col]])) +
          geom_point(size = 1.6, alpha = 0.7) +
          stat_ellipse(aes(group = .data[[group_col]]), level = 0.95,
                       linewidth = 0.6, show.legend = FALSE) +
          scale_colour_manual(values = group_palette, drop = FALSE) +
          labs(title = paste(lab, "PCoA"),
               x = paste0("PCoA1 (", round(pc$prop[1] * 100, 2), "%)"),
               y = paste0("PCoA2 (", round(pc$prop[2] * 100, 2), "%)"),
               colour = group_col) + theme_paper()
        save_plot(g, out("02_beta/Figures", paste0(lab, "_PCoA.pdf")), 7, 6)
      }
    }

    # NMDS
    nd <- read_nmds(metric)
    if (!is.null(nd)) {
      nd <- inner_join(nd, md, by = "sample-id")
      if (nrow(nd) >= 3) {
        g <- ggplot(nd, aes(NMDS1, NMDS2, colour = .data[[group_col]])) +
          geom_point(size = 1.6, alpha = 0.7) +
          stat_ellipse(aes(group = .data[[group_col]]), level = 0.95,
                       linewidth = 0.6, show.legend = FALSE) +
          scale_colour_manual(values = group_palette, drop = FALSE) +
          labs(title = paste(lab, "NMDS"), colour = group_col) + theme_paper()
        save_plot(g, out("02_beta/Figures", paste0(lab, "_NMDS.pdf")), 7, 6)
      }
    } else if (recalculate_nmds) {
      nm <- tryCatch(metaMDS(d, k = 2, trymax = 100, autotransform = FALSE, trace = 0),
                     error = function(e) NULL)
      if (!is.null(nm)) {
        nd <- as.data.frame(scores(nm, display = "sites"))
        nd$`sample-id` <- rownames(nd)
        nd <- inner_join(nd, md, by = "sample-id")
        g <- ggplot(nd, aes(NMDS1, NMDS2, colour = .data[[group_col]])) +
          geom_point(size = 1.6, alpha = 0.7) +
          scale_colour_manual(values = group_palette, drop = FALSE) +
          labs(title = paste(lab, "NMDS; stress", round(nm$stress, 3)),
               colour = group_col) + theme_paper()
        save_plot(g, out("02_beta/Figures", paste0(lab, "_NMDS_recalculated.pdf")), 7, 6)
      }
    }

    # PERMANOVA: group + Region（Region 缺失时只检验 group）
    terms <- unique(c(group_col, if (region_col %in% names(md)) region_col))
    md_model <- md[complete.cases(md[, terms, drop = FALSE]), , drop = FALSE]
    if (nrow(md_model) >= 3 && nlevels(droplevels(as.factor(md_model[[group_col]]))) >= 2) {
      dm <- as.dist(dmat[md_model[["sample-id"]], md_model[["sample-id"]], drop = FALSE])
      form <- as.formula(paste("dm ~", paste(sprintf("`%s`", terms), collapse = " + ")))
      ad <- tryCatch(adonis2(form, data = md_model, permutations = permutations, by = "margin"),
                     error = function(e) { warning("PERMANOVA ", lab, ": ", conditionMessage(e)); NULL })
      if (!is.null(ad)) {
        z <- as.data.frame(ad); z$Term <- rownames(z); z$Distance <- lab
        permanova_results[[lab]] <- z
      }
    }

    # PERMDISP: 组内离散度
    disp_group <- droplevels(as.factor(md[[group_col]]))
    if (nlevels(disp_group) >= 2) {
      bd <- tryCatch(betadisper(d, disp_group, type = "median"), error = function(e) NULL)
      if (!is.null(bd)) {
        bt <- tryCatch(permutest(bd, permutations = permutations), error = function(e) NULL)
        if (!is.null(bt)) {
          z <- as.data.frame(bt$tab); z$Term <- rownames(z); z$Distance <- lab
          permdisp_results[[lab]] <- z
        }
        dd <- data.frame(sample_id = ids, Group = disp_group, DistanceToCentroid = bd$distances)
        write_table(dd, out("02_beta/Tables", paste0(lab, "_PERMDISP_distances.tsv")))
        g <- ggplot(dd, aes(Group, DistanceToCentroid, fill = Group)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.75) +
          geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
          scale_fill_manual(values = group_palette, drop = FALSE) +
          labs(title = paste(lab, "PERMDISP"), x = group_col,
               y = "Distance to group centroid") +
          theme_paper() + theme(legend.position = "none")
        save_plot(g, out("02_beta/Figures", paste0(lab, "_PERMDISP.pdf")), 7, 5.5)
      }
    }
  }

  if (length(permanova_results)) {
    z <- bind_rows(permanova_results)
    pc <- grep("^Pr\\(>F\\)$", names(z), value = TRUE)
    if (length(pc) == 1) z <- z %>% group_by(Term) %>% mutate(FDR = p.adjust(.data[[pc]], "BH")) %>% ungroup()
    write_table(z, out("02_beta/Tables/PERMANOVA_all_distances.tsv"))
  }
  if (length(permdisp_results)) {
    z <- bind_rows(permdisp_results)
    pc <- grep("^Pr\\(>F\\)$", names(z), value = TRUE)
    if (length(pc) == 1) z$FDR <- p.adjust(z[[pc]], "BH")
    write_table(z, out("02_beta/Tables/PERMDISP_all_distances.tsv"))
  }
  sf <- p("05_Diversity", "nmds", "NMDS_stress_summary.tsv")
  if (file.exists(sf)) write_table(read_tsv(sf), out("02_beta/Tables/NMDS_stress_summary.tsv"))
}

############################ 7. COMPOSITION ####################################

if (run_composition) {
  plot_composition <- function(mat, rank, top_n) {
    rel <- relative_abundance(mat)
    av <- sort(colMeans(rel), decreasing = TRUE)
    top <- head(names(av)[names(av) != "Unclassified"], top_n)
    keep <- intersect(top, colnames(rel))
    dat <- rel[, keep, drop = FALSE]
    other <- rowSums(rel[, setdiff(colnames(rel), keep), drop = FALSE])
    if (any(other > 0)) dat <- cbind(dat, Other = other)
    long <- as.data.frame(dat) %>%
      rownames_to_column("sample-id") %>%
      pivot_longer(-`sample-id`, names_to = "Taxon", values_to = "Abundance") %>%
      left_join(meta_rare, by = "sample-id") %>%
      filter(!is.na(.data[[group_col]]))
    summ <- long %>%
      group_by(.data[[group_col]], Taxon) %>%
      summarise(MeanAbundance = mean(Abundance), .groups = "drop")
    ord <- summ %>% group_by(Taxon) %>% summarise(m = mean(MeanAbundance), .groups = "drop") %>%
      arrange(m) %>% pull(Taxon)
    summ$Taxon <- factor(summ$Taxon, levels = ord)
    g <- ggplot(summ, aes(.data[[group_col]], MeanAbundance, fill = Taxon)) +
      geom_col(width = 0.78) +
      scale_y_continuous(labels = percent_format()) +
      labs(x = group_col, y = "Mean relative abundance",
           title = paste(rank, "community composition"), fill = rank) +
      theme_paper() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
    save_plot(g, out("03_composition", paste0(rank, "_composition.pdf")), 9, 6)
    write_table(summ, out("11_tables", paste0(rank, "_composition_by_", group_col, ".tsv")))
  }
  plot_composition(phylum_mat, "Phylum", top_phylum_n)
  plot_composition(genus_mat, "Genus", top_genus_n)
}

############################ 8. GENUS DIFFERENTIAL: KW + BH ####################

if (run_genus_differential) {
  sample_ids <- intersect(rownames(genus_mat), meta_rare[["sample-id"]])
  if (length(sample_ids) < 3) stop("属水平丰度表与元数据匹配的样本数不足")
  genus_mat2 <- genus_mat[sample_ids, , drop = FALSE]
  meta_genus <- meta_rare[match(sample_ids, meta_rare[["sample-id"]]), , drop = FALSE]

  genus_rel <- sweep(genus_mat2, 1, rowSums(genus_mat2), "/")
  genus_rel[!is.finite(genus_rel)] <- 0

  long <- as.data.frame(genus_rel, check.names = FALSE) %>%
    rownames_to_column("sample-id") %>%
    pivot_longer(cols = -`sample-id`, names_to = "Genus", values_to = "Relative_abundance") %>%
    left_join(meta_genus %>% select(`sample-id`, all_of(group_col)), by = "sample-id") %>%
    filter(!is.na(.data[[group_col]]))

  stats <- long %>%
    group_by(Genus) %>%
    summarise(
      prevalence = mean(Relative_abundance > 0, na.rm = TRUE),
      mean_abundance = mean(Relative_abundance, na.rm = TRUE),
      p_value = kw_p(Relative_abundance, .data[[group_col]]),
      .groups = "drop"
    ) %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR, p_value)

  write_table(stats, out("04_genus_differential/Genus_KW_BH_all.tsv"))

  top <- stats %>% filter(!is.na(FDR), FDR < fdr_cutoff) %>%
    slice_head(n = top_diff_n) %>% pull(Genus)

  if (length(top) > 0) {
    hd <- long %>%
      filter(Genus %in% top) %>%
      group_by(.data[[group_col]], Genus) %>%
      summarise(mean_abundance = mean(Relative_abundance, na.rm = TRUE), .groups = "drop")
    hd$Genus <- factor(hd$Genus, levels = rev(top))
    g <- ggplot(hd, aes(.data[[group_col]], Genus, fill = mean_abundance)) +
      geom_tile(colour = "white") +
      scale_fill_gradient(low = "white", high = "#B2182B",
                          labels = percent_format(accuracy = 0.01)) +
      labs(x = group_col, y = "Genus", fill = "Mean relative abundance",
           title = paste0("Top ", length(top), " genera (KW-BH ranking)")) +
      theme_paper()
    save_plot(g, out("04_genus_differential/Top_genus_heatmap.pdf"),
              8, max(6, min(12, 0.25 * length(top) + 3)))
  } else {
    message("没有属达到当前 FDR 阈值：", fdr_cutoff, "，跳过 Top genus 热图")
  }

  volc <- stats %>%
    mutate(neglog10FDR = -log10(pmax(FDR, .Machine$double.xmin)),
           Significance = ifelse(!is.na(FDR) & FDR < fdr_cutoff, "FDR < 0.05", "Not significant"))
  g <- ggplot(volc, aes(log10(mean_abundance + 1e-8), neglog10FDR, colour = Significance)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_colour_manual(values = c("FDR < 0.05" = "#D73027", "Not significant" = "grey65")) +
    labs(x = "log10(mean relative abundance)", y = "-log10(BH-FDR)",
         title = "Genus-level exploratory differential analysis") +
    theme_paper()
  save_plot(g, out("04_genus_differential/Genus_volcano.pdf"), 7, 6)
}

############################ 9. CORE MICROBIOME ################################

if (run_core) {
  rel <- relative_abundance(genus_mat)
  md <- meta_rare[match(rownames(rel), meta_rare[["sample-id"]]), , drop = FALSE]
  core_by_group <- lapply(group_levels, function(gr) {
    ids <- which(as.character(md[[group_col]]) == gr)
    if (!length(ids)) return(NULL)
    data.frame(Genus = colnames(rel), Group = gr,
               Prevalence = colMeans(rel[ids, , drop = FALSE] > 0),
               Mean_relative_abundance = colMeans(rel[ids, , drop = FALSE]))
  })
  core_df <- bind_rows(core_by_group)
  write_table(core_df, out("05_core/Core_prevalence_by_group.tsv"))

  core_genus <- core_df %>%
    group_by(Genus) %>%
    summarise(Max_group_prevalence = max(Prevalence),
              Mean_relative_abundance = mean(Mean_relative_abundance), .groups = "drop") %>%
    filter(Max_group_prevalence >= core_prevalence) %>%
    arrange(desc(Mean_relative_abundance)) %>%
    slice_head(n = core_top_n)
  write_table(core_genus, out("05_core/Core_genera_threshold.tsv"))

  if (nrow(core_genus)) {
    hp <- core_df %>% filter(Genus %in% core_genus$Genus) %>%
      mutate(Genus = factor(Genus, levels = rev(core_genus$Genus)))
    g <- ggplot(hp, aes(Group, Genus, fill = Prevalence)) +
      geom_tile(colour = "white") +
      scale_fill_viridis_c(limits = c(0, 1), labels = percent_format()) +
      labs(x = group_col, y = "Genus", fill = "Prevalence",
           title = paste0("Core genera: any-group prevalence \u2265 ", core_prevalence * 100, "%")) +
      theme_paper()
    save_plot(g, out("05_core/Core_prevalence_heatmap.pdf"),
              8, max(6, min(14, 0.25 * nrow(core_genus) + 3)))

    ha <- core_df %>% filter(Genus %in% core_genus$Genus) %>%
      mutate(Genus = factor(Genus, levels = rev(core_genus$Genus)))
    g2 <- ggplot(ha, aes(Group, Genus, fill = Mean_relative_abundance)) +
      geom_tile(colour = "white") +
      scale_fill_viridis_c(labels = percent_format(accuracy = 0.01)) +
      labs(x = group_col, y = "Genus", fill = "Mean relative abundance") +
      theme_paper()
    save_plot(g2, out("05_core/Core_abundance_heatmap.pdf"),
              8, max(6, min(14, 0.25 * nrow(core_genus) + 3)))
  }

  sens <- bind_rows(lapply(c(0.5, 0.6, 0.7, 0.8, 0.9), function(th) {
    core_df %>% group_by(Genus) %>% summarise(max_prev = max(Prevalence), .groups = "drop") %>%
      summarise(threshold = th, n_core = sum(max_prev >= th))
  }))
  write_table(sens, out("05_core/Core_threshold_sensitivity.tsv"))
}

############################ 10. VENN ##########################################

if (run_venn) {
  if (!requireNamespace("ggvenn", quietly = TRUE)) {
    warning("Venn 图需要 ggvenn，已跳过。")
  } else {
    make_venn <- function(mat, rank) {
      sets <- lapply(group_levels, function(gr) {
        ids <- meta_rare[["sample-id"]][as.character(meta_rare[[group_col]]) == gr]
        ids <- intersect(ids, rownames(mat))
        if (!length(ids)) return(character(0))
        present <- colnames(mat)[colSums(mat[ids, , drop = FALSE] > 0) > 0]
        present[!grepl("unclassified|unknown|uncultured|^NA$", present, ignore.case = TRUE)]
      })
      names(sets) <- group_levels
      write_table(data.frame(Group = names(sets), Taxa_count = lengths(sets)),
                  out("06_venn", paste0(rank, "_Venn_group_counts.tsv")))

      if (length(sets) >= 2 && length(sets) <= 4 && all(lengths(sets) > 0)) {
        plt <- ggvenn::ggvenn(
          sets,
          fill_color = unname(group_palette[names(sets)]),
          fill_alpha = 0.35, stroke_color = "black", stroke_size = 0.7,
          set_name_color = "black", set_name_size = 4,
          text_color = "black", text_size = 3.5,
          show_percentage = FALSE, show_elements = FALSE
        ) +
          labs(title = paste0("Shared and unique ", rank, " among ", group_col, " groups")) +
          theme_void() +
          theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13, colour = "black"),
                plot.margin = margin(12, 12, 12, 12))
        save_plot(plt, out("06_venn", paste0(rank, "_Venn.pdf")), 8, 7)

        all_taxa <- unique(unlist(sets))
        membership <- data.frame(Taxon = all_taxa)
        for (gr in names(sets)) membership[[gr]] <- membership$Taxon %in% sets[[gr]]
        membership$Combination <- apply(membership[, names(sets), drop = FALSE], 1,
                                        function(z) paste(names(sets)[z], collapse = ";"))
        write_table(membership, out("06_venn", paste0(rank, "_Venn_membership.tsv")))
      } else if (length(sets) > 4) {
        warning(rank, " Venn 图超过 4 组，ggvenn 不绘图；组计数已输出。")
      } else {
        warning(rank, " 存在空分组，无法绘制 Venn 图；组计数已输出。")
      }
    }
    make_venn(phylum_mat, "Phylum")
    make_venn(genus_mat, "Genus")
  }
}

############################ 11. LEfSe ########################################

if (run_lefse) {
  # 自包含的 LEfSe 风格分析（Kruskal-Wallis + 组间 LDA 效应量），
  # 不依赖 microbiomeMarker / phyloseq，避免其空 marker_table 的 validObject 报错。
  mat <- as.matrix(genus_mat); storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0; mat[mat < 0] <- 0
  md <- meta_rare[match(rownames(mat), meta_rare[["sample-id"]]), , drop = FALSE]
  ok <- !is.na(md[[group_col]]) & !is.na(md[["sample-id"]])
  mat <- mat[ok, , drop = FALSE]; md <- md[ok, , drop = FALSE]
  md[[group_col]] <- droplevels(as.factor(md[[group_col]]))
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  mat <- mat[, !grepl("^(Unclassified|Other|NA)$", colnames(mat), ignore.case = TRUE), drop = FALSE]

  zero_lib <- rowSums(mat) <= 0
  if (any(zero_lib)) {
    mat <- mat[!zero_lib, , drop = FALSE]
    md <- md[!zero_lib, , drop = FALSE]
    md[[group_col]] <- droplevels(md[[group_col]])
  }

  if (nrow(mat) < 3 || ncol(mat) < 2 || nlevels(md[[group_col]]) < 2) {
    warning("LEfSe 有效样本/属/分组不足，跳过。")
  } else {
    grp <- md[[group_col]]
    classes <- levels(grp)

    # 1) KW + BH（全属）
    rel <- sweep(mat, 1, rowSums(mat), "/")
    kw_long <- as.data.frame(rel) %>%
      rownames_to_column("sample-id") %>%
      pivot_longer(-`sample-id`, names_to = "Genus", values_to = "RelAbundance") %>%
      left_join(md %>% select(`sample-id`, all_of(group_col)), by = "sample-id")
    lefse_kw <- kw_long %>%
      group_by(Genus) %>%
      summarise(n = n(), P_value = kw_p(RelAbundance, .data[[group_col]]),
                .groups = "drop") %>%
      mutate(FDR_BH = p.adjust(P_value, "BH")) %>%
      arrange(FDR_BH, P_value)
    write_table(lefse_kw, out("07_LEfSe/LEfSe_Genus_KW_BH_all.tsv"))

    # 2) KW 显著属：富集组 + LDA 效应量
    cpm <- sweep(mat, 1, rowSums(mat), "/") * 1e6   # CPM 尺度，与 LEfSe 的 LDA 效应量尺度一致
    sig_gen <- lefse_kw$Genus[!is.na(lefse_kw$P_value) & lefse_kw$P_value <= lefse_kw_cutoff]

    if (length(sig_gen) >= 2) {
      # 富集组 = 该属 CPM 均值最高的组
      grp_mean <- t(vapply(sig_gen, function(g) {
        vapply(classes, function(c) mean(cpm[grp == c, g]), numeric(1))
      }, numeric(length(classes))))
      colnames(grp_mean) <- classes
      enrich_grp <- classes[apply(grp_mean, 1, which.max)]

      # LDA 效应量 = log10(1 + 最高组 CPM 均值 - 其余组 CPM 均值)，始终为正
      gap <- vapply(seq_along(sig_gen), function(i) {
        top  <- grp_mean[i, enrich_grp[i]]
        rest <- mean(grp_mean[i, classes != enrich_grp[i]])
        top - rest
      }, numeric(1))
      ef_lda <- log10(1 + pmax(gap, 0))

      markers <- data.frame(
        feature = sig_gen,
        enrich_group = enrich_grp,
        ef_lda = as.numeric(ef_lda),
        stringsAsFactors = FALSE
      ) %>%
        left_join(lefse_kw %>% select(Genus, P_value, FDR_BH),
                  by = c("feature" = "Genus")) %>%
        rename(pvalue = P_value, padj = FDR_BH) %>%
        arrange(enrich_group, desc(ef_lda))
    } else {
      markers <- data.frame(feature = character(), enrich_group = character(),
                            ef_lda = numeric(), pvalue = numeric(), padj = numeric())
    }
    write_table(markers, out("07_LEfSe/LEfSe_Genus_all_markers.tsv"))

    # 3) 最终差异属：BH-FDR < fdr_cutoff 且 |LDA| > lefse_final_lda
    lefse_final <- markers %>%
      filter(!is.na(padj), padj < fdr_cutoff, is.finite(ef_lda),
             abs(ef_lda) > lefse_final_lda) %>%
      arrange(padj, desc(abs(ef_lda)))
    write_table(lefse_final, out("07_LEfSe/LEfSe_Genus_BH0.05_LDA4.tsv"))

    # 4) LDA 柱状图
    if (nrow(lefse_final) > 0) {
      pd <- lefse_final %>%
        mutate(Enriched_in = enrich_group, LDA_magnitude = abs(ef_lda),
               Genus = reorder(feature, LDA_magnitude))
      g <- ggplot(pd, aes(LDA_magnitude, Genus, fill = Enriched_in)) +
        geom_col(width = 0.72) +
        scale_fill_manual(values = group_palette, drop = FALSE) +
        labs(x = "Absolute LDA score", y = "Genus", fill = "Enriched in",
             title = paste("LEfSe-style", group_col)) +
        theme_paper()
      save_plot(g, out("07_LEfSe/LEfSe_Genus_LDA_barplot.pdf"),
                8, max(5, 0.28 * nrow(pd) + 2))
    }

    # 5) 摘要
    write_table(data.frame(
      parameter = c("group_col", "taxonomic_level", "method", "normalization",
                    "KW_cutoff", "LDA_final_cutoff", "significant_genera", "n_samples_used"),
      value = c(group_col, "Genus", "LEfSe-style (KW + group-mean LDA score)", "CPM",
                lefse_kw_cutoff, lefse_final_lda, nrow(lefse_final), nrow(mat))
    ), out("07_LEfSe/LEfSe_analysis_summary.tsv"))
  }
}

############################ 12. FUNCTIONAL PROFILES ###########################

# KEGG 通路 koID -> 通路名映射。
# Tax4Fun2 导出的 pathway 表行名只有 koXXXXX（通路编号）而没有通路名称，
# 故这里内置常用 KEGG 通路的名称，用于给结果表和热图加上可读名称。
kegg_pathway_names <- c(
  "ko00010" = "Glycolysis / Gluconeogenesis",
  "ko00020" = "Citrate cycle (TCA cycle)",
  "ko00030" = "Pentose phosphate pathway",
  "ko00040" = "Pentose and glucuronate interconversions",
  "ko00051" = "Fructose and mannose metabolism",
  "ko00052" = "Galactose metabolism",
  "ko00053" = "Ascorbate and aldarate metabolism",
  "ko00061" = "Fatty acid biosynthesis",
  "ko00062" = "Fatty acid elongation",
  "ko00071" = "Fatty acid degradation",
  "ko00072" = "Synthesis and degradation of ketone bodies",
  "ko00100" = "Steroid biosynthesis",
  "ko00120" = "Primary bile acid biosynthesis",
  "ko00121" = "Secondary bile acid biosynthesis",
  "ko00130" = "Ubiquinone and other terpenoid-quinone biosynthesis",
  "ko00190" = "Oxidative phosphorylation",
  "ko00195" = "Photosynthesis",
  "ko00196" = "Photosynthesis - antenna proteins",
  "ko00220" = "Arginine biosynthesis",
  "ko00230" = "Purine metabolism",
  "ko00240" = "Pyrimidine metabolism",
  "ko00250" = "Alanine, aspartate and glutamate metabolism",
  "ko00260" = "Glycine, serine and threonine metabolism",
  "ko00261" = "Monobactam biosynthesis",
  "ko00270" = "Cysteine and methionine metabolism",
  "ko00280" = "Valine, leucine and isoleucine degradation",
  "ko00290" = "Valine, leucine and isoleucine biosynthesis",
  "ko00300" = "Lysine biosynthesis",
  "ko00310" = "Lysine degradation",
  "ko00330" = "Arginine and proline metabolism",
  "ko00340" = "Histidine metabolism",
  "ko00350" = "Tyrosine metabolism",
  "ko00360" = "Phenylalanine metabolism",
  "ko00362" = "Benzoate degradation",
  "ko00380" = "Tryptophan metabolism",
  "ko00400" = "Phenylalanine, tyrosine and tryptophan biosynthesis",
  "ko00410" = "beta-Alanine metabolism",
  "ko00430" = "Taurine and hypotaurine metabolism",
  "ko00440" = "Phosphonate and phosphinate metabolism",
  "ko00450" = "Selenocompound metabolism",
  "ko00460" = "Cyanoamino acid metabolism",
  "ko00471" = "D-Glutamine and D-glutamate metabolism",
  "ko00472" = "D-Arginine and D-ornithine metabolism",
  "ko00473" = "D-Alanine metabolism",
  "ko00480" = "Glutathione metabolism",
  "ko00500" = "Starch and sucrose metabolism",
  "ko00510" = "N-Glycan biosynthesis",
  "ko00511" = "Other glycan degradation",
  "ko00520" = "Amino sugar and nucleotide sugar metabolism",
  "ko00521" = "Streptomycin biosynthesis",
  "ko00540" = "Lipopolysaccharide biosynthesis",
  "ko00550" = "Peptidoglycan biosynthesis",
  "ko00561" = "Glycerolipid metabolism",
  "ko00562" = "Inositol phosphate metabolism",
  "ko00564" = "Glycerophospholipid metabolism",
  "ko00565" = "Ether lipid metabolism",
  "ko00590" = "Arachidonic acid metabolism",
  "ko00591" = "Linoleic acid metabolism",
  "ko00592" = "alpha-Linolenic acid metabolism",
  "ko00600" = "Sphingolipid metabolism",
  "ko00620" = "Pyruvate metabolism",
  "ko00630" = "Glyoxylate and dicarboxylate metabolism",
  "ko00640" = "Propanoate metabolism",
  "ko00650" = "Butanoate metabolism",
  "ko00660" = "C5-Branched dibasic acid metabolism",
  "ko00670" = "One carbon pool by folate",
  "ko00680" = "Methane metabolism",
  "ko00710" = "Carbon fixation in photosynthetic organisms",
  "ko00720" = "Carbon fixation pathways in prokaryotes",
  "ko00730" = "Thiamine metabolism",
  "ko00740" = "Riboflavin metabolism",
  "ko00750" = "Vitamin B6 metabolism",
  "ko00760" = "Nicotinate and nicotinamide metabolism",
  "ko00770" = "Pantothenate and CoA biosynthesis",
  "ko00780" = "Biotin metabolism",
  "ko00785" = "Lipoic acid metabolism",
  "ko00790" = "Folate biosynthesis",
  "ko00860" = "Porphyrin metabolism",
  "ko00900" = "Terpenoid backbone biosynthesis",
  "ko00906" = "Carotenoid biosynthesis",
  "ko00910" = "Nitrogen metabolism",
  "ko00920" = "Sulfur metabolism",
  "ko00970" = "Aminoacyl-tRNA biosynthesis",
  "ko00980" = "Metabolism of xenobiotics by cytochrome P450",
  "ko01040" = "Biosynthesis of unsaturated fatty acids",
  "ko01053" = "Biosynthesis of siderophore group nonribosomal peptides",
  "ko01100" = "Metabolic pathways",
  "ko01110" = "Biosynthesis of secondary metabolites",
  "ko01120" = "Microbial metabolism in diverse environments",
  "ko01130" = "Biosynthesis of antibiotics",
  "ko01200" = "Carbon metabolism",
  "ko01210" = "2-Oxocarboxylic acid metabolism",
  "ko01212" = "Fatty acid metabolism",
  "ko01220" = "Degradation of aromatic compounds",
  "ko01230" = "Biosynthesis of amino acids",
  "ko01240" = "Biosynthesis of cofactors",
  "ko01250" = "Biosynthesis of nucleotide sugars",
  "ko01260" = "Carbon metabolism in prokaryotes",
  "ko02010" = "ABC transporters",
  "ko02020" = "Two-component system",
  "ko02024" = "Quorum sensing",
  "ko02025" = "Biofilm formation - Vibrio cholerae",
  "ko02026" = "Biofilm formation - Pseudomonas aeruginosa",
  "ko02030" = "Bacterial chemotaxis",
  "ko02040" = "Flagellar assembly",
  "ko02060" = "Phosphotransferase system (PTS)",
  "ko03008" = "Ribosome biogenesis in eukaryotes",
  "ko03010" = "Ribosome",
  "ko03013" = "RNA transport",
  "ko03015" = "mRNA surveillance pathway",
  "ko03018" = "RNA degradation",
  "ko03020" = "RNA polymerase",
  "ko03022" = "Basal transcription factors",
  "ko03030" = "DNA replication",
  "ko03040" = "Spliceosome",
  "ko03050" = "Proteasome",
  "ko03060" = "Protein export",
  "ko03070" = "Bacterial secretion system",
  "ko03410" = "Base excision repair",
  "ko03420" = "Nucleotide excision repair",
  "ko03430" = "Mismatch repair",
  "ko03440" = "Homologous recombination",
  "ko03450" = "Non-homologous end-joining",
  "ko03460" = "Fanconi anemia pathway",
  "ko04011" = "MAPK signaling pathway - yeast",
  "ko04112" = "Cell cycle - Caulobacter",
  "ko04122" = "Sulfur relay system",
  "ko04146" = "Peroxisome",
  "ko04626" = "Plant-pathogen interaction",
  "ko01501" = "beta-Lactam resistance",
  "ko01502" = "Vancomycin resistance",
  "ko01503" = "Cationic antimicrobial peptide (CAMP) resistance",
  "ko02061" = "Protein export (bacteria)",
  "ko05110" = "Vibrio cholerae infection",
  "ko05111" = "Biofilm formation - Vibrio",
  "ko05120" = "Epithelial cell signaling in Helicobacter pylori infection",
  "ko00140" = "Steroid hormone biosynthesis",
  "ko00232" = "Caffeine metabolism",
  "ko00253" = "Tetracycline biosynthesis",
  "ko00281" = "Geraniol degradation",
  "ko00311" = "Penicillin and cephalosporin biosynthesis",
  "ko00331" = "Clavulanic acid biosynthesis",
  "ko00332" = "Carbapenem biosynthesis",
  "ko00333" = "Prodigiosin biosynthesis",
  "ko00361" = "Chlorocyclohexane and chlorobenzene degradation",
  "ko00363" = "Bisphenol degradation",
  "ko00364" = "Fluorobenzoate degradation",
  "ko00365" = "Furfural degradation",
  "ko00401" = "Novobiocin biosynthesis",
  "ko00404" = "Staurosporine biosynthesis",
  "ko00405" = "Phenazine biosynthesis",
  "ko00513" = "Various types of N-glycan biosynthesis",
  "ko00514" = "Other types of O-glycan biosynthesis",
  "ko00515" = "Mannose type O-glycan biosynthesis",
  "ko00522" = "Biosynthesis of 12-, 14- and 16-membered macrolides",
  "ko00523" = "Polyketide sugar unit biosynthesis",
  "ko00524" = "Neomycin, kanamycin and gentamicin biosynthesis",
  "ko00525" = "Acarbose and validamycin biosynthesis",
  "ko00531" = "Glycosaminoglycan degradation",
  "ko00532" = "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",
  "ko00534" = "Glycosaminoglycan biosynthesis - heparan sulfate / heparin",
  "ko00563" = "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
  "ko00571" = "Lipoarabinomannan (LAM) biosynthesis",
  "ko00601" = "Glycosphingolipid biosynthesis - lacto and neolacto series"
)

# 获取完整的 KEGG 通路名称（koID -> 通路名）。
# 优先联网从 KEGG 官方 REST API 抓取完整列表（覆盖全部通路，跟分组变量无关）；
# 离线或抓取失败时，回退到上面的内置常用映射。
get_kegg_pathway_names <- function() {
  builtin <- kegg_pathway_names
  old_timeout <- getOption("timeout")
  options(timeout = 40)
  on.exit(options(timeout = old_timeout), add = TRUE)
  x <- tryCatch(
    readLines("https://rest.kegg.jp/list/pathway/ko", warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || length(x) < 10) {
    message("未能联网获取 KEGG 完整通路名称，使用内置映射（", length(builtin), " 条）。")
    return(builtin)
  }
  parts <- strsplit(x, "\t", fixed = TRUE)
  ids <- vapply(parts, function(z) z[1], character(1))
  nms <- vapply(parts, function(z) if (length(z) >= 2) z[2] else "", character(1))
  out <- nms; names(out) <- ids
  out <- out[nzchar(ids) & nzchar(nms)]
  out <- out[!duplicated(names(out))]   # 去重，防止同名导致向量索引异常
  if (length(out) >= length(builtin)) {
    message("已联网获取 KEGG 完整通路名称（", length(out), " 条）。")
    out
  } else {
    builtin
  }
}

# 统一采用 代码3 的协变量校正方法：
#   lm(log10(相对丰度% + pseudocount) ~ group + Region) vs lm(~ Region)
#   嵌套 anova 做全局检验；emmeans 给校正均值；pairs 做两两比较；BH-FDR 校正。
read_function_table <- function(file, id_candidates, meta_samples) {
  x <- read.delim(file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
                  quote = "", comment.char = "", row.names = NULL)
  if (!nzchar(names(x)[1]) || is.na(names(x)[1])) names(x)[1] <- "Feature_ID"
  id_col <- intersect(id_candidates, names(x))[1]
  if (is.na(id_col)) id_col <- names(x)[1]
  sample_cols <- intersect(names(x), meta_samples)
  if (length(sample_cols) < 3) stop("功能表匹配到的样本列少于 3 个：", file)
  ids <- make.unique(as.character(x[[id_col]]))
  mm <- as.matrix(x[, sample_cols, drop = FALSE])
  storage.mode(mm) <- "numeric"; mm[!is.finite(mm)] <- 0
  rownames(mm) <- ids
  t(mm)   # samples x features
}

analyse_function_table <- function(mat, label, top_n = top_function_n, pairwise_mode = "all", name_map = NULL) {
  # 给功能项 ID 附加可读名称（例如通路 koID -> 通路名），无映射时保留 ID。
  # 写成向量化版本，供 vapply（标量）和 mutate（整列向量）两种情况共用。
  fun_label <- function(f) {
    f <- as.character(f)
    if (is.null(name_map)) return(f)
    lbl <- name_map[f]               # 按名查找，未匹配为 NA
    out <- f
    ok <- !is.na(lbl)
    out[ok] <- paste0(lbl[ok], " (", f[ok], ")")
    out
  }
  ids <- intersect(meta_full[["sample-id"]], rownames(mat))
  if (length(ids) < 3) { warning(label, " 匹配样本不足。"); return(NULL) }
  md <- meta_full[match(ids, meta_full[["sample-id"]]), , drop = FALSE]
  mm <- mat[ids, , drop = FALSE]
  mm <- mm[, colSums(mm, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(mm) == 0) { warning(label, " 没有非零功能项，跳过。"); return(NULL) }

  md[[group_col]] <- droplevels(as.factor(md[[group_col]]))
  has_region <- region_col %in% names(md) && nlevels(droplevels(as.factor(md[[region_col]]))) >= 1
  if (has_region) md[[region_col]] <- droplevels(as.factor(md[[region_col]]))

  # 相对丰度 (%) + log10
  rel <- sweep(mm, 1, rowSums(mm), "/") * 100
  rel[!is.finite(rel)] <- 0

  full_terms <- unique(c(group_col, if (has_region) region_col))
  full_formula <- as.formula(paste0("y ~ ", paste(sprintf("`%s`", full_terms), collapse = " + ")))
  red_formula  <- if (has_region) {
    as.formula(paste0("y ~ `", region_col, "`"))
  } else {
    y ~ 1
  }

  features <- colnames(rel)
  global_results <- list(); failed_results <- list()

  # 设计矩阵秩（与响应无关，只计算一次）
  design_dat <- data.frame(y = rep(0, nrow(md)))
  design_dat[[group_col]] <- md[[group_col]]
  if (has_region) design_dat[[region_col]] <- md[[region_col]]
  full_rank <- qr(model.matrix(full_formula, data = design_dat))$rank

  # ---- Pass 1: 全局检验 lm(full) vs lm(reduced) + anova ----
  for (i in seq_along(features)) {
    f <- features[i]
    y <- log10(rel[, f] + pseudocount)
    dat <- data.frame(y = y)
    dat[[group_col]] <- md[[group_col]]
    if (has_region) dat[[region_col]] <- md[[region_col]]

    full <- tryCatch(lm(full_formula, data = dat, na.action = na.fail), error = function(e) NULL)
    red  <- tryCatch(lm(red_formula,  data = dat, na.action = na.fail), error = function(e) NULL)
    if (is.null(full) || is.null(red)) { failed_results[[f]] <- "模型拟合失败"; next }
    if (full$rank < full_rank || df.residual(full) <= 0) { failed_results[[f]] <- "共线性/自由度不足"; next }

    a <- tryCatch(anova(red, full), error = function(e) NULL)
    if (is.null(a) || nrow(a) < 2) { failed_results[[f]] <- "嵌套比较失败"; next }
    global_results[[f]] <- data.frame(Function = f, F_value = a$F[2],
                                      P_value = a$`Pr(>F)`[2])
    if (i %% 500 == 0) cat("  [", label, "] 全局检验 ", i, "/", length(features), "\n", sep = "")
  }

  global_df <- bind_rows(global_results)
  if (nrow(global_df) == 0) { warning(label, " 无有效全局检验结果。"); return(NULL) }
  global_df$FDR <- p.adjust(global_df$P_value, method = "BH")
  global_df <- global_df %>% arrange(FDR)

  # ---- Pass 2: emmeans 校正均值 + 两两比较（只对需要的功能项）----
  fun_mean <- colMeans(rel)
  top_funs <- names(sort(fun_mean, decreasing = TRUE))[seq_len(min(top_n, length(fun_mean)))]
  sig_feat <- global_df$Function[!is.na(global_df$FDR) & global_df$FDR < fdr_cutoff]
  emm_features <- if (pairwise_mode == "all") global_df$Function else unique(c(top_funs, sig_feat))

  adjusted_results <- list(); pairwise_results <- list()
  for (i in seq_along(emm_features)) {
    f <- emm_features[i]
    y <- log10(rel[, f] + pseudocount)
    dat <- data.frame(y = y)
    dat[[group_col]] <- md[[group_col]]
    if (has_region) dat[[region_col]] <- md[[region_col]]
    full <- tryCatch(lm(full_formula, data = dat, na.action = na.fail), error = function(e) NULL)
    if (is.null(full)) next
    emm <- tryCatch(emmeans::emmeans(full, specs = as.formula(paste0("~ `", group_col, "`"))),
                    error = function(e) NULL)
    if (is.null(emm)) next
    emm_df <- as.data.frame(emm)
    adjusted_results[[f]] <- data.frame(
      Function = f, Group = as.character(emm_df[[group_col]]),
      Adjusted_mean_log10 = emm_df$emmean, SE = emm_df$SE,
      Lower_CL = emm_df$lower.CL, Upper_CL = emm_df$upper.CL
    )
    if (pairwise_mode == "all" || f %in% sig_feat) {
      pw <- tryCatch(pairs(emm, adjust = "none"), error = function(e) NULL)
      if (!is.null(pw)) {
        pw_df <- as.data.frame(pw)
        pairwise_results[[f]] <- data.frame(
          Function = f, Contrast = pw_df$contrast, Estimate = pw_df$estimate,
          SE = pw_df$SE, P_value = pw_df$p.value
        )
      }
    }
  }

  safe <- gsub("[^A-Za-z0-9_]+", "_", label)
  global_df$Name <- vapply(global_df$Function, fun_label, character(1))
  write_table(global_df, out("08_function", paste0(safe, "_global_test_BH.tsv")))

  adj_df <- bind_rows(adjusted_results)
  if (nrow(adj_df)) {
    adj_df$Name <- vapply(adj_df$Function, fun_label, character(1))
    write_table(adj_df, out("08_function", paste0(safe, "_adjusted_means.tsv")))
  }

  pw_df <- bind_rows(pairwise_results)
  if (nrow(pw_df)) {
    pw_df$FDR <- p.adjust(pw_df$P_value, method = "BH")
    pw_df <- pw_df %>% arrange(FDR)
    pw_df$Name <- vapply(pw_df$Function, fun_label, character(1))
    write_table(pw_df, out("08_function", paste0(safe, "_pairwise_BH.tsv")))
  }

  failed_df <- data.frame(Function = names(failed_results),
                          Reason = unlist(failed_results, use.names = FALSE))
  if (nrow(failed_df)) write_table(failed_df, out("08_function", paste0(safe, "_model_failed.tsv")))

  # ---- 热图矩阵（校正均值 -> 行 z-score）----
  make_heatmap_matrix <- function(funs) {
    d <- adj_df %>% filter(Function %in% funs) %>%
      select(Function, Group, Adjusted_mean_log10)
    if (!nrow(d)) return(NULL)
    wide <- d %>% pivot_wider(names_from = Group, values_from = Adjusted_mean_log10)
    m <- as.data.frame(wide); rownames(m) <- m$Function; m$Function <- NULL
    m <- as.matrix(m); storage.mode(m) <- "numeric"
    keep <- apply(m, 1, function(z) all(is.finite(z)) && sd(z) > 0)
    m <- m[keep, , drop = FALSE]
    if (!nrow(m)) return(NULL)
    z <- t(scale(t(m))); z[!is.finite(z)] <- 0
    z
  }
  plot_function_heatmap <- function(zmat, title_text, filename) {
    if (is.null(zmat) || nrow(zmat) == 0) return(NULL)
    pd <- as.data.frame(zmat) %>%
      rownames_to_column("Function") %>%
      pivot_longer(-Function, names_to = "Group", values_to = "Zscore") %>%
      left_join(global_df %>% select(Function, FDR), by = "Function") %>%
      mutate(Function_label = fun_label(Function),
             Function_label = ifelse(!is.na(FDR) & FDR < fdr_cutoff,
                                     paste0(Function_label, " *"), Function_label))
    pd$Function_label <- factor(pd$Function_label, levels = rev(unique(pd$Function_label)))
    g <- ggplot(pd, aes(Group, Function_label, fill = Zscore)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "Row Z-score") +
      labs(title = title_text, x = group_col, y = NULL) +
      theme_paper() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 7))
    save_plot(g, out("08_function", filename), 9, max(6, min(14, 0.28 * nrow(zmat) + 3)))
    g
  }

  plot_function_heatmap(make_heatmap_matrix(top_funs),
                        paste0(label, ": top ", length(top_funs), " functions"),
                        paste0(safe, "_Top_functions_heatmap.pdf"))

  sig_funs <- global_df %>% filter(!is.na(FDR), FDR < fdr_cutoff) %>% pull(Function)
  if (length(sig_funs) > 0) {
    plot_function_heatmap(make_heatmap_matrix(sig_funs),
                          paste0(label, ": significant functions (BH-FDR < ", fdr_cutoff, ")"),
                          paste0(safe, "_Significant_functions_heatmap.pdf"))
  }
  invisible(global_df)
}

if (run_functional_analysis) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    warning("功能分析需要 emmeans，已跳过。")
  } else {
    suppressPackageStartupMessages(library(emmeans))
    fap     <- p("06_Function/FAPROTAX/functional_table.tsv")
    ko      <- p("06_Function/Tax4Fun2/Tax4Fun2_KO.tsv")
    pathway <- p("06_Function/Tax4Fun2/Tax4Fun2_pathway.tsv")
    for (f in c(fap, ko, pathway)) if (!file.exists(f)) warning("功能预测文件不存在：", f)

    if (file.exists(fap))
      analyse_function_table(read_function_table(fap, c("group", "Function"), meta_full[["sample-id"]]),
                             "FAPROTAX", 15, pairwise_mode = "all")
    if (file.exists(pathway)) {
      pathway_names <- get_kegg_pathway_names()   # 完整映射（联网）或内置回退
      analyse_function_table(read_function_table(pathway, c("pathway", "Pathway", "Description", "ko"), meta_full[["sample-id"]]),
                             "Tax4Fun2_Pathway", 15, pairwise_mode = "all",
                             name_map = pathway_names)
    }
    if (file.exists(ko))
      analyse_function_table(read_function_table(ko, c("KO", "ko", "Feature"), meta_full[["sample-id"]]),
                             "Tax4Fun2_KO", 15, pairwise_mode = "significant")
  }
}

############################ 13. CONTINUOUS TRAITS #############################

if (run_continuous_traits) {
  traits <- intersect(c("TL", "Max_length"), names(meta_rare))
  cont <- list()
  for (lab in names(distance_names)) {
    dm <- tryCatch(read_distance(distance_names[[lab]]), error = function(e) NULL)
    if (is.null(dm)) next
    for (v in traits) {
      md <- meta_rare[is.finite(meta_rare[[v]]), , drop = FALSE]
      if (region_col %in% names(md)) md <- md[!is.na(md[[region_col]]), , drop = FALSE]
      ids <- intersect(md[["sample-id"]], rownames(dm))
      if (length(ids) < 10 || length(unique(md[[v]])) < 3) next
      md <- md[match(ids, md[["sample-id"]]), , drop = FALSE]
      d <- as.dist(dm[ids, ids, drop = FALSE])
      terms <- unique(c(v, if (region_col %in% names(md)) region_col))
      ff <- as.formula(paste("d ~", paste(sprintf("`%s`", terms), collapse = " + ")))
      rr <- tryCatch(adonis2(ff, data = md, permutations = permutations, by = "margin"),
                     error = function(e) NULL)
      if (!is.null(rr)) {
        z <- as.data.frame(rr); z$Term <- rownames(z); z$Trait <- v; z$Distance <- lab
        cont[[paste(lab, v)]] <- z
      }
    }
  }
  if (length(cont)) write_table(bind_rows(cont), out("09_continuous_traits/continuous_trait_PERMANOVA.tsv"))

  corrs <- list()
  if (run_alpha && nrow(alpha_long) > 0 && length(traits) > 0) {
    for (v in traits) for (m in unique(alpha_long$metric)) {
      z <- alpha_long %>% filter(is.finite(.data[[v]]), metric == m)
      if (nrow(z) >= 5 && length(unique(z[[v]])) >= 3) {
        ct <- suppressWarnings(cor.test(z$value, z[[v]], method = "spearman", exact = FALSE))
        corrs[[paste(m, v)]] <- data.frame(Metric = m, Trait = v, n = nrow(z),
                                           rho = unname(ct$estimate), P_value = ct$p.value)
      }
    }
  }
  if (length(corrs)) write_table(bind_rows(corrs) %>% mutate(FDR = p.adjust(P_value, "BH")),
                                 out("09_continuous_traits/alpha_trait_Spearman.tsv"))
}

############################ 14. ANCOM-BC2 ####################################

if (run_ancombc2) {
  if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
    warning("ANCOMBC 未安装，跳过 ANCOM-BC2。可运行 BiocManager::install('ANCOMBC') 安装。")
  } else {
    suppressPackageStartupMessages(library(ANCOMBC))

    raw_biom_file <- p("02_Abundance/BIOM/otu_table_no_eukaryota.biom")
    tax_file <- p("03_Taxonomy_Sequences/taxonomy_no_eukaryota.tsv")

    # 1) 直接读 HDF5，得到 feature x sample 稀疏矩阵
    raw0 <- read_hdf5_biom_sparse(raw_biom_file)
    raw_asv <- Matrix::t(raw0)   # sample x feature

    meta_an <- read_tsv(p("01_Metadata/metadata_full.tsv")) %>% distinct(`sample-id`, .keep_all = TRUE)
    tax_an <- read_tsv(tax_file)
    fc <- intersect(c("Feature ID", "Feature.ID", "feature-id", "FeatureID", "ASV"), names(tax_an))[1]
    if (is.na(fc)) fc <- names(tax_an)[1]
    tc <- intersect(c("Taxon", "taxonomy", "Taxonomy"), names(tax_an))[1]
    if (is.na(tc)) stop("taxonomy 缺少 Taxon 列。")
    tx <- tax_an %>%
      transmute(FeatureID = as.character(.data[[fc]]), Taxon = as.character(.data[[tc]])) %>%
      distinct(FeatureID, .keep_all = TRUE)
    tx$Genus <- vapply(tx$Taxon, function(x) {
      z <- split_tax(x); h <- z[grepl("^g__", z, ignore.case = TRUE)]
      a <- if (length(h)) sub("^g__", "", tail(h, 1), ignore.case = TRUE) else NA_character_
      if (is.na(a) || !nzchar(a) ||
          grepl("unclassified|uncultured|unknown|metagenome|incertae_sedis", a, ignore.case = TRUE))
        NA_character_ else trimws(a)
    }, character(1))

    # 用 FeatureID 名查找属名，避免按位置索引导致与 BIOM 列顺序不一致时引入 NA
    genus_lookup <- setNames(tx$Genus, tx$FeatureID)
    genus_vec <- unname(genus_lookup[colnames(raw_asv)])   # 每个 ASV 的属名，未匹配为 NA
    keep <- which(!is.na(genus_vec))
    if (!length(keep)) stop("原始 BIOM 与 genus taxonomy 没有匹配 ASV。")
    gen <- genus_vec[keep]
    gen_levels <- sort(unique(gen))
    gmap <- Matrix::sparseMatrix(i = match(gen, gen_levels), j = seq_along(keep), x = 1,
                                 dims = c(length(gen_levels), length(keep)))
    genus_counts <- raw_asv[, keep, drop = FALSE] %*% Matrix::t(gmap)  # sample x genus
    rownames(genus_counts) <- rownames(raw_asv)
    colnames(genus_counts) <- gen_levels

    meta_an <- meta_an[match(rownames(genus_counts), meta_an[["sample-id"]]), , drop = FALSE]
    good <- !is.na(meta_an[[group_col]]) & !is.na(meta_an[["sample-id"]])
    if (region_col %in% names(meta_an)) good <- good & !is.na(meta_an[[region_col]])
    meta_an <- meta_an[good, , drop = FALSE]
    genus_counts <- genus_counts[good, , drop = FALSE]

    lib <- Matrix::rowSums(genus_counts)
    keep_s <- lib >= 1000
    genus_counts <- genus_counts[keep_s, , drop = FALSE]
    meta_an <- meta_an[keep_s, , drop = FALSE]
    if (nrow(genus_counts) < 10) stop("ANCOM-BC2 过滤测序深度后样本数不足（要求 genus 总计数>=1000）。")

    genus_counts <- genus_counts[, Matrix::colSums(genus_counts) > 0, drop = FALSE]

    meta_an[[group_col]] <- droplevels(as.factor(meta_an[[group_col]]))
    if (!is.null(reference_level) && reference_level %in% levels(meta_an[[group_col]]))
      meta_an[[group_col]] <- relevel(meta_an[[group_col]], ref = reference_level)
    if (nlevels(meta_an[[group_col]]) < 2) stop("ANCOM-BC2 分组变量不足两个水平。")
    if (region_col %in% names(meta_an)) meta_an[[region_col]] <- droplevels(as.factor(meta_an[[region_col]]))

    formula_terms <- unique(c(group_col, if (region_col %in% names(meta_an)) region_col))
    # ANCOM-BC2 的 data_sanity_check 按字面比较公式中的变量名，不能带反引号
    fix_formula <- paste(formula_terms, collapse = " + ")
    mm_design <- model.matrix(as.formula(paste("~", fix_formula)), data = meta_an)
    if (qr(mm_design)$rank < ncol(mm_design))
      stop("ANCOM-BC2 模型矩阵秩亏；请检查分组与 Region 是否完全混杂或存在空水平。")

    counts_for_ancom <- t(as.matrix(genus_counts))
    storage.mode(counts_for_ancom) <- "numeric"
    rownames(meta_an) <- meta_an[["sample-id"]]
    meta_an[["sample-id"]] <- NULL   # ANCOM-BC2 用行名作为样本 ID，去掉冗余列
    if (!identical(colnames(counts_for_ancom), rownames(meta_an)))
      stop("ANCOM-BC2 count matrix 与 metadata 样本顺序不一致。")

    set.seed(123)
    ancom_result <- ANCOMBC::ancombc2(
      data = counts_for_ancom, taxa_are_rows = TRUE, meta_data = meta_an,
      fix_formula = fix_formula, rand_formula = NULL, p_adj_method = "BH",
      pseudo = 0, pseudo_sens = TRUE, prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
      group = group_col, struc_zero = TRUE, neg_lb = TRUE, alpha = fdr_cutoff,
      n_cl = 1, verbose = TRUE, global = TRUE, pairwise = TRUE,
      dunnet = FALSE, trend = FALSE,
      iter_control = list(tol = 1e-5, max_iter = 100, verbose = FALSE),
      em_control = list(tol = 1e-5, max_iter = 100),
      lme_control = NULL,
      mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
      trend_control = NULL
    )
    saveRDS(ancom_result, out("10_ANCOMBC2/ANCOMBC2_full_result.rds"))

    for (nm in c("res", "res_global", "res_pair", "zero_ind", "pseudo_sens_tab"))
      if (!is.null(ancom_result[[nm]]))
        write_table(as.data.frame(ancom_result[[nm]]), out("10_ANCOMBC2", paste0("ANCOMBC2_", nm, ".tsv")))

    # 全局显著
    if (!is.null(ancom_result$res_global)) {
      gdf <- as.data.frame(ancom_result$res_global)
      if (all(c("q_val", "diff_abn") %in% names(gdf))) {
        gsig <- gdf[!is.na(gdf$q_val) & gdf$q_val < fdr_cutoff & gdf$diff_abn, , drop = FALSE]
        write_table(gsig, out("10_ANCOMBC2/ANCOMBC2_global_significant.tsv"))
      }
    }

    # 主效应 LFC（group 系数，动态识别列名）
    if (!is.null(ancom_result$res)) {
      res <- as.data.frame(ancom_result$res)
      taxon_col <- "taxon"
      lfc_cols <- grep("^lfc_", names(res), value = TRUE)
      lfc_cols <- lfc_cols[!grepl("(Intercept)", lfc_cols, fixed = TRUE)]
      if (region_col %in% names(meta_an) && region_col != group_col)
        lfc_cols <- lfc_cols[!startsWith(lfc_cols, paste0("lfc_", region_col))]

      if (length(lfc_cols)) {
        lfc_long <- do.call(rbind, lapply(lfc_cols, function(lc) {
          suffix <- sub("^lfc_", "", lc)
          data.frame(taxon = res[[taxon_col]], comparison = suffix,
                     lfc = res[[lc]],
                     se = if (paste0("se_", suffix) %in% names(res)) res[[paste0("se_", suffix)]] else NA_real_,
                     q = if (paste0("q_", suffix) %in% names(res)) res[[paste0("q_", suffix)]] else NA_real_,
                     diff = if (paste0("diff_", suffix) %in% names(res)) res[[paste0("diff_", suffix)]] else NA,
                     stringsAsFactors = FALSE)
        }))
        write_table(lfc_long, out("10_ANCOMBC2/ANCOMBC2_Group_LFC_all.tsv"))
        lfc_sig <- lfc_long[!is.na(lfc_long$q) & lfc_long$q < fdr_cutoff &
                              !is.na(lfc_long$lfc) & abs(lfc_long$lfc) >= 1 &
                              lfc_long$diff %in% TRUE, , drop = FALSE]
        write_table(lfc_sig, out("10_ANCOMBC2/ANCOMBC2_Group_LFC_significant.tsv"))
        if (nrow(lfc_sig)) {
          lfc_top <- do.call(rbind, lapply(unique(lfc_sig$comparison), function(cmp) {
            x <- lfc_sig[lfc_sig$comparison == cmp, , drop = FALSE]
            head(x[order(x$q, -abs(x$lfc)), , drop = FALSE], 15)
          }))
          write_table(lfc_top, out("10_ANCOMBC2/ANCOMBC2_Group_LFC_Top15.tsv"))
          p_lfc <- ggplot(lfc_top, aes(lfc, reorder(taxon, lfc))) +
            geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
            geom_point(size = 2.6, colour = "#2878B5") +
            facet_wrap(~comparison, scales = "free_y") +
            labs(title = "ANCOM-BC2: top group-associated genera",
                 x = "Log fold change", y = "Genus") +
            theme_paper() + theme(axis.text.y = element_text(size = 9))
          save_plot(p_lfc, out("10_ANCOMBC2/ANCOMBC2_Group_LFC_Top15.pdf"), 10, 7)
        }
      }
    }

    # 两两比较
    if (!is.null(ancom_result$res_pair)) {
      pr <- as.data.frame(ancom_result$res_pair)
      lfc_cols <- grep("^lfc_", names(pr), value = TRUE)
      lfc_cols <- lfc_cols[!grepl("(Intercept)", lfc_cols, fixed = TRUE)]
      if (region_col %in% names(meta_an) && region_col != group_col)
        lfc_cols <- lfc_cols[!startsWith(lfc_cols, paste0("lfc_", region_col))]
      if (length(lfc_cols)) {
        pw_long <- do.call(rbind, lapply(lfc_cols, function(lc) {
          suffix <- sub("^lfc_", "", lc)
          data.frame(taxon = pr$taxon, comparison = suffix,
                     lfc = pr[[lc]],
                     q = if (paste0("q_", suffix) %in% names(pr)) pr[[paste0("q_", suffix)]] else NA_real_,
                     diff = if (paste0("diff_", suffix) %in% names(pr)) pr[[paste0("diff_", suffix)]] else NA,
                     diff_robust = if (paste0("diff_robust_", suffix) %in% names(pr)) pr[[paste0("diff_robust_", suffix)]] else NA,
                     stringsAsFactors = FALSE)
        }))
        write_table(pw_long, out("10_ANCOMBC2/ANCOMBC2_Pairwise_all.tsv"))
        pw_sig <- pw_long[!is.na(pw_long$q) & pw_long$q < fdr_cutoff &
                            pw_long$diff_robust %in% TRUE, , drop = FALSE]
        write_table(pw_sig, out("10_ANCOMBC2/ANCOMBC2_Pairwise_robust_significant.tsv"))
        if (nrow(pw_sig)) {
          pw_top <- do.call(rbind, lapply(unique(pw_sig$comparison), function(cmp) {
            x <- pw_sig[pw_sig$comparison == cmp, , drop = FALSE]
            head(x[order(x$q, -abs(x$lfc)), , drop = FALSE], 15)
          }))
          write_table(pw_top, out("10_ANCOMBC2/ANCOMBC2_Pairwise_Top15.tsv"))
          pw_top$rowid <- paste(pw_top$taxon, pw_top$comparison, sep = "___")
          p_pw <- ggplot(pw_top, aes(lfc, reorder(rowid, lfc))) +
            geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
            geom_point(size = 2.4, colour = "#2878B5") +
            facet_wrap(~comparison, scales = "free_y", ncol = 2) +
            scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
            labs(title = "ANCOM-BC2: pairwise group comparisons",
                 x = "Log fold change", y = "Genus") +
            theme_paper() + theme(axis.text.y = element_text(size = 8))
          save_plot(p_pw, out("10_ANCOMBC2/ANCOMBC2_Pairwise_Top15.pdf"), 12, 10)
        }
      }
    }

    write_table(data.frame(
      parameter = c("group_col", "reference_level", "fixed_formula",
                    "prevalence_cutoff", "minimum_genus_library_size",
                    "p_adjustment", "pairwise_adjustment", "struc_zero",
                    "pseudo_sensitivity", "ANCOMBC_version", "n_samples", "n_genera"),
      value = c(group_col, levels(meta_an[[group_col]])[1], fix_formula, 0.10, 1000,
                "BH", "Holm", "TRUE", "TRUE", as.character(packageVersion("ANCOMBC")),
                nrow(meta_an), ncol(counts_for_ancom))
    ), out("10_ANCOMBC2/ANCOMBC2_settings.tsv"))
  }
}

############################ 15. FINAL REPORT ##################################

settings <- data.frame(
  parameter = c("package_dir", "group_col", "region_col", "reference_level",
                "permutations", "rarefied_samples", "rarefied_ASVs",
                "genus_number", "phylum_number", "functional_metadata_samples"),
  value = c(package_dir, group_col, region_col,
            ifelse(is.null(reference_level), "automatic", reference_level),
            permutations, nrow(meta_rare), ncol(otu), ncol(genus_mat),
            ncol(phylum_mat), nrow(meta_full))
)
write_table(settings, out("12_QC/analysis_settings.tsv"))

files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
write_table(data.frame(file = substring(files, nchar(out_dir) + 2),
                       size_bytes = file.info(files)$size),
            out("12_QC/results_manifest.tsv"))

cat("\n====================================================\n")
cat("Fish 16S meta-analysis 完成\n")
cat("分组变量：", group_col, "\n")
cat("Region 协变量：", region_col, "\n")
cat("稀释样本数：", nrow(meta_rare), "\n")
cat("稀释 ASV 数：", ncol(otu), "\n")
cat("结果目录：", out_dir, "\n")
cat("注意：KW 属差异分析为探索性检验；LEfSe 为独立分析；\n")
cat("ANCOM-BC2 使用未稀释原始计数并按 Genus 聚合；\n")
cat("功能分析为协变量校正后的 lm + emmeans；\n")
cat("PERMDISP 检验组内离散度，不等同于 PERMANOVA。\n")
cat("====================================================\n")
