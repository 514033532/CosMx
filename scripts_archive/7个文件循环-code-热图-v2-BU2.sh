library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# -------------------------- 1. 通用配置（所有文件共用） --------------------------
# 1.1 定义样本顺序和分组信息
sample_order <- c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6")

# 分组方案1：DN vs NC
group_info1 <- data.frame(
  Sample = sample_order,
  Group = c("DN", "DN", "DN", "DN", "DN", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = FALSE,
  row.names = sample_order
)

# 分组方案2：DN1/DN2 vs NC
group_info2 <- data.frame(
  Sample = sample_order,
  Group = c("DN2", "DN2", "DN1", "DN2", "DN1", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = FALSE,
  row.names = sample_order
)

# 1.2 定义颜色方案
heat_colors <- colorRampPalette(c(
  "#006400", "#228B22", "#32CD32", "#90EE90",  # 下调：深绿→浅绿
  "#F5F5F5",                                     # 中性色
  "#FFB6C1", "#FF6347", "#DC143C", "#B22222"    # 上调：浅红→深红
))(100)
annotation_colors1 <- list(Group = c("DN" = "#E41A1C", "NC" = "#377EB8"))
annotation_colors2 <- list(Group = c("DN1" = "#E41A1C", "DN2" = "#FF7F00", "NC" = "#377EB8"))
change_colors <- list(Change = c("up in case" = "#E41A1C", "down in case" = "#377EB8"))

# 1.3 修正文件路径：拆分「差异列表文件」和「热图矩阵文件」（格式正确+分类清晰）
# --- 差异列表文件（含Change/Significance列，_list后缀）---
list_file_paths <- list(
  GeneMoran_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/spatial_distr/moran/Group2_DN-VS-NC/GeneMoran_DE_list_in_Group2_DN-VS-NC.csv",
  CellType_prop_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/proportion/Group2_DN-VS-NC/CellType_prop_DE_list_in_Group2_DN-VS-NC.csv"
)

# --- 热图矩阵文件（无Change/Significance列，_heatmap后缀）---
heatmap_file_paths <- list(
  GeneMoran_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/spatial_distr/moran/Group2_DN-VS-NC/GeneMoran_DE_heatmap_in_Group2_DN-VS-NC.csv",
  CellType_prop_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/proportion/Group2_DN-VS-NC/CellType_prop_DE_heatmap_in_Group2_DN-VS-NC.csv"
)

# -------------------------- 2. 差异列表文件处理函数（含Change/Significance） --------------------------
plot_heatmap_from_list <- function(file_path, prefix) {
  cat(paste("开始处理差异列表文件：", file_path, "\n"))
  data <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  # 检查必需列
  required_cols <- c("Change", "Significance")
  if (!all(required_cols %in% colnames(data))) {
    missing <- setdiff(required_cols, colnames(data))
    stop(paste(prefix, "文件缺少必需列：", paste(missing, collapse = ", ")))
  }
  
  # 筛选显著差异数据
  sig_data <- data %>% filter(Significance == "Sig")
  if (nrow(sig_data) == 0) {
    cat(paste(prefix, "文件无显著差异数据，跳过绘图\n"))
    return()
  }
  
  # 筛选前20个上调/下调
  up_data <- sig_data %>% filter(Change == "up in case") %>% head(20)
  down_data <- sig_data %>% filter(Change == "down in case") %>% head(20)
  plot_data <- rbind(up_data, down_data)
  
  if (nrow(plot_data) < 1) {
    cat(paste(prefix, "文件无符合条件的上调/下调数据，跳过绘图\n"))
    return()
  }
  
  # 提取样本列并标准化
  sample_data <- plot_data[, sample_order, drop = FALSE]
  matrix_data <- as.matrix(sample_data)
  matrix_scaled <- t(scale(t(matrix_data)))  # 行标准化
  
  # 行注释（上调/下调）
  annotation_row <- data.frame(
    Change = factor(plot_data$Change, levels = c("up in case", "down in case"))
  )
  rownames(annotation_row) <- rownames(plot_data)
  
  # 输出目录
  output_dir <- dirname(file_path)
  setwd(output_dir)
  
  # 绘制热图（分组1：DN vs NC）
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_row = annotation_row,
    annotation_col = group_info1["Group"],  # 移除无效的drop=FALSE
    annotation_colors = c(annotation_colors1, change_colors),
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Expression (DN vs NC)"),
    filename = paste0(prefix, "_heatmap_group1.pdf"),
    width = 14,
    height = max(10, nrow(plot_data) * 0.5),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14
  )
  
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_row = annotation_row,
    annotation_col = group_info1["Group"],
    annotation_colors = c(annotation_colors1, change_colors),
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Expression (DN vs NC)"),
    filename = paste0(prefix, "_heatmap_group1.png"),
    width = 14,
    height = max(10, nrow(plot_data) * 0.5),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14,
    dpi = 300
  )
  
  # 绘制热图（分组2：DN1/DN2 vs NC）
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_row = annotation_row,
    annotation_col = group_info2["Group"],  # 移除无效的drop=FALSE
    annotation_colors = c(annotation_colors2, change_colors),
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Expression (DN1, DN2 vs NC)"),
    filename = paste0(prefix, "_heatmap_group2.pdf"),
    width = 14,
    height = max(10, nrow(plot_data) * 0.5),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14
  )
  
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_row = annotation_row,
    annotation_col = group_info2["Group"],
    annotation_colors = c(annotation_colors2, change_colors),
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Expression (DN1, DN2 vs NC)"),
    filename = paste0(prefix, "_heatmap_group2.png"),
    width = 14,
    height = max(10, nrow(plot_data) * 0.5),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14,
    dpi = 300
  )
  
  cat(paste(prefix, "文件热图已保存到：", output_dir, "\n\n"))
}

# -------------------------- 3. 热图矩阵文件处理函数（无Change/Significance） --------------------------
plot_heatmap_from_matrix <- function(file_path, prefix) {
  cat(paste("开始处理热图矩阵文件：", file_path, "\n"))
  data <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  # 仅保留样本列（按sample_order筛选，避免列名不匹配）
  sample_cols <- intersect(sample_order, colnames(data))
  if (length(sample_cols) == 0) {
    stop(paste(prefix, "文件无匹配的样本列！"))
  }
  plot_data <- data[, sample_cols, drop = FALSE]
  
  if (nrow(plot_data) < 1) {
    cat(paste(prefix, "文件无数据，跳过绘图\n"))
    return()
  }
  
  # 行标准化
  matrix_data <- as.matrix(plot_data)
  matrix_scaled <- t(scale(t(matrix_data)))
  
  # 输出目录
  output_dir <- dirname(file_path)
  setwd(output_dir)
  
  # 绘制热图（分组1：DN vs NC，无行注释）
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_col = group_info1["Group"],
    annotation_colors = annotation_colors1,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Expression (DN vs NC)"),
    filename = paste0(prefix, "_heatmap_group1.pdf"),
    width = 14,
    height = max(10, nrow(plot_data) * 0.5),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14
  )
  
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_col = group_info1["Group"],
    annotation_colors = annotation_colors1,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Expression (DN vs NC)"),
    filename = paste0(prefix, "_heatmap_group1.png"),
    width = 14,
    height = max(10, nrow(plot_data) * 0.5),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14,
    dpi = 300
  )
  
  # 绘制热图（分组2：DN1/DN2 vs NC，无行注释）
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_col = group_info2["Group"],
    annotation_colors = annotation_colors2,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Expression (DN1, DN2 vs NC)"),
    filename = paste0(prefix, "_heatmap_group2.pdf"),
    width = 14,
    height = max(10, nrow(plot_data) * 0.5),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14
  )
  
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_col = group_info2["Group"],
    annotation_colors = annotation_colors2,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Expression (DN1, DN2 vs NC)"),
    filename = paste0(prefix, "_heatmap_group2.png"),
    width = 14,
    height = max(10, nrow(plot_data) * 0.5),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14,
    dpi = 300
  )
  
  cat(paste(prefix, "文件热图已保存到：", output_dir, "\n\n"))
}

# -------------------------- 4. 批量执行 --------------------------
# 处理差异列表文件（_list后缀）
if (length(list_file_paths) > 0) {
  for (i in seq_along(list_file_paths)) {
    plot_heatmap_from_list(list_file_paths[[i]], names(list_file_paths)[i])
  }
}

# 处理热图矩阵文件（_heatmap后缀）
if (length(heatmap_file_paths) > 0) {
  for (i in seq_along(heatmap_file_paths)) {
    plot_heatmap_from_matrix(heatmap_file_paths[[i]], names(heatmap_file_paths)[i])
  }
}

cat("所有文件处理完成！\n")
