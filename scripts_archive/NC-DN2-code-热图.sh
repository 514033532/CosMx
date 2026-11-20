library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# 1. 设置工作目录
setwd("/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_colocate/celltype_enrich/Group1_DN2-VS-NC/")

# 2. 读取数据
file_path <- "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_colocate/celltype_enrich/Group1_DN2-VS-NC/CellEnrich_DE_list_in_Group1_DN2-VS-NC.csv"
data <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# 3. 定义样本信息：仅保留 DN2 和 NC 组（核心筛选）
sample_info <- data.frame(
  Sample = c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6"),
  Group = c("DN2", "DN2", "DN1", "DN2", "DN1", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = FALSE
)
# 筛选 DN2 和 NC 样本，避免 DN1 干扰
sample_info_filtered <- sample_info %>% filter(Group %in% c("DN2", "NC"))
sample_order_filtered <- sample_info_filtered$Sample
group_info2 <- sample_info_filtered
rownames(group_info2) <- group_info2$Sample  # 行名设为样本名，匹配热图列

# 4. 检查必需列
required_cols <- c("Change", "Significance")
if (!all(required_cols %in% colnames(data))) {
  missing <- setdiff(required_cols, colnames(data))
  stop(paste("数据缺少必需列：", paste(missing, collapse = ", ")))
}

# 5. 筛选显著差异 LR 对
sig_data <- data %>% filter(Significance == "Sig")
if (nrow(sig_data) == 0) {
  stop("无显著差异（Significance = 'Sig'）的 LR 对，无法绘图")
}

# 6. 筛选 up/down LR 对（各取前20）
up_data <- sig_data %>% filter(Change == "up in case") %>% head(20)
down_data <- sig_data %>% filter(Change == "down in case") %>% head(20)
plot_data <- rbind(up_data, down_data)
if (nrow(plot_data) < 1) {
  stop("无符合条件的 up/down LR 对，无法绘图")
}

# 7. 提取样本数据（仅 DN2/NC）
missing_samples <- setdiff(sample_order_filtered, colnames(plot_data))
if (length(missing_samples) > 0) {
  warning(paste("缺少样本列，已跳过：", paste(missing_samples, collapse = ", ")))
  sample_order_filtered <- intersect(sample_order_filtered, colnames(plot_data))
  if (length(sample_order_filtered) == 0) {
    stop("无有效样本列，终止运行")
  }
}
sample_data <- plot_data[, sample_order_filtered, drop = FALSE]

# 8. 转换矩阵并标准化（处理非数值数据）
matrix_data <- as.matrix(sample_data)
matrix_data <- apply(matrix_data, 2, function(x) as.numeric(as.character(x)))  # 确保数值型
matrix_scaled <- t(scale(t(matrix_data)))  # 行标准化
rownames(matrix_scaled) <- rownames(plot_data)

# 9. 定义注释信息
annotation_row <- data.frame(
  Change = factor(plot_data$Change, levels = c("up in case", "down in case")),
  row.names = rownames(plot_data)
)
annotation_col <- group_info2[, "Group", drop = FALSE]  # 列注释：DN2/NC 分组

# 10. 修正调色板：用 RdBu 基础调色板 + rev() 实现反转（解决原错误）
# 先获取 RdBu 调色板的11个基础颜色，再用 rev() 反转（低表达蓝→高表达红）
heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)  

# 11. 定义注释颜色（与之前风格统一）
annotation_colors <- list(
  Group = c("DN2" = "#FF7F00", "NC" = "#006400"),  # DN2=橙，NC=深绿
  Change = c("up in case" = "#E41A1C", "down in case" = "#377EB8")  # up=红，down=蓝
)

# 12. 绘制 PDF 热图（发表用）
pheatmap(
  mat = matrix_scaled,
  scale = "none",
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 15,
  treeheight_col = 15,
  fontsize = 12,
  fontsize_row = 10,
  fontsize_col = 12,
  color = heat_colors,
  main = "LR Pairs Expression (DN2 vs NC)",
  filename = "LR_pair_heatmap_DN2_vs_NC.pdf",
  width = 12,
  height = ifelse(nrow(plot_data) > 20, 14, 10),
  border_color = NA,
  cluster_cols = TRUE,
  cellwidth = 20,
  cellheight = 12
)

# 13. 绘制 PNG 热图（展示用）
pheatmap(
  mat = matrix_scaled,
  scale = "none",
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 15,
  treeheight_col = 15,
  fontsize = 12,
  fontsize_row = 10,
  fontsize_col = 12,
  color = heat_colors,
  main = "LR Pairs Expression (DN2 vs NC)",
  filename = "LR_pair_heatmap_DN2_vs_NC.png",
  width = 12,
  height = ifelse(nrow(plot_data) > 20, 14, 10),
  border_color = NA,
  cluster_cols = TRUE,
  cellwidth = 20,
  cellheight = 12,
  dpi = 300
)

cat("热图已生成！文件：LR_pair_heatmap_DN2_vs_NC.pdf / .png\n")
