library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# 1. 设置工作目录（确保路径正确）
setwd("/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_colocate/celltype_enrich/Group1_DN2-VS-DN1/")

# 2. 读取数据（注意：原文件路径是 DN2-VS-NC，需确认是否与分析目标一致）
file_path <- "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_colocate/celltype_enrich/Group1_DN2-VS-NC/CellEnrich_DE_list_in_Group1_DN2-VS-NC.csv"
data <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# 3. 定义样本信息：仅保留 DN2 和 NC 组，确保“样本名-分组”一一对应（核心修正）
# 先明确所有样本的原始分组，再筛选 DN2 和 NC
sample_info <- data.frame(
  Sample = c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6"),
  Group = c("DN2", "DN2", "DN1", "DN2", "DN1", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = FALSE
)

# 筛选：仅保留 DN2 和 NC 组的样本
sample_info_filtered <- sample_info %>% filter(Group %in% c("DN2", "NC"))
sample_order_filtered <- sample_info_filtered$Sample  # 筛选后的样本顺序
group_info2 <- sample_info_filtered  # 最终用于注释的分组信息（仅 DN2/NC）
rownames(group_info2) <- group_info2$Sample  # 设置行名为样本名，匹配热图列名

# 4. 检查必需列是否存在（避免报错）
required_cols <- c("Change", "Significance")
if (!all(required_cols %in% colnames(data))) {
  missing <- setdiff(required_cols, colnames(data))
  stop(paste("数据中缺少必需列：", paste(missing, collapse = ", ")))
}

# 5. 筛选显著差异的 LR 对（Significance = "Sig"）
sig_data <- data %>% filter(Significance == "Sig")
if (nrow(sig_data) == 0) {
  stop("没有显著差异（Significance = 'Sig'）的 LR 对，无法绘制热图")
}

# 6. 分别筛选“up in case”和“down in case”的 LR 对（各取前20个，不足则取全部）
up_data <- sig_data %>% filter(Change == "up in case") %>% head(20)
down_data <- sig_data %>% filter(Change == "down in case") %>% head(20)

# 合并数据（确保至少有1行数据）
plot_data <- rbind(up_data, down_data)
if (nrow(plot_data) < 1) {
  stop("没有符合条件的 LR 对（up/down in case），无法绘制热图")
}

# 7. 提取筛选后样本的表达数据（仅保留 DN2/NC 样本）
# 检查数据中是否存在筛选后的样本列
missing_samples <- setdiff(sample_order_filtered, colnames(plot_data))
if (length(missing_samples) > 0) {
  warning(paste("数据中缺少以下样本列，已跳过：", paste(missing_samples, collapse = ", ")))
  sample_order_filtered <- intersect(sample_order_filtered, colnames(plot_data))
  if (length(sample_order_filtered) == 0) {
    stop("没有找到有效的样本列，无法继续")
  }
}
sample_data <- plot_data[, sample_order_filtered, drop = FALSE]

# 8. 转换为矩阵并进行行标准化（使每行 LR 对的表达可比较）
matrix_data <- as.matrix(sample_data)
# 处理可能的非数值型数据（避免标准化报错）
matrix_data <- apply(matrix_data, 2, function(x) as.numeric(as.character(x)))
matrix_scaled <- t(scale(t(matrix_data)))  # 行标准化（按 LR 对标准化）
rownames(matrix_scaled) <- rownames(plot_data)  # 恢复 LR 对名称

# 9. 定义注释信息（行：LR 对的变化类型；列：样本分组）
# 行注释：标记 LR 对是“up”还是“down”
annotation_row <- data.frame(
  Change = factor(plot_data$Change, levels = c("up in case", "down in case")),
  row.names = rownames(plot_data)
)

# 列注释：样本的 DN2/NC 分组（已在步骤3筛选）
annotation_col <- group_info2[, "Group", drop = FALSE]  # 仅保留 Group 列

# 10. 定义注释颜色（统一配色，符合发表格式）
annotation_colors <- list(
  Group = c("DN2" = "#FF7F00", "NC" = "#006400"),  # DN2=橙色，NC=深绿色
  Change = c("up in case" = "#E41A1C", "down in case" = "#377EB8")  # up=红色，down=蓝色
)

# 11. 定义热图颜色梯度（蓝→白→红，标准化后的值对应）
heat_colors <- colorRampPalette(brewer.pal(11, "RdBu_r"))(100)  # RdBu_r 是红→蓝反转，即低表达蓝、高表达红

# 12. 绘制热图（PDF格式，适合发表）
pheatmap(
  mat = matrix_scaled,
  scale = "none",  # 已提前行标准化，此处无需重复
  annotation_row = annotation_row,  # 行注释（LR 变化类型）
  annotation_col = annotation_col,  # 列注释（样本分组：DN2/NC）
  annotation_colors = annotation_colors,  # 注释颜色
  show_rownames = TRUE,  # 显示 LR 对名称
  show_colnames = TRUE,  # 显示样本名
  treeheight_row = 15,  # 行聚类树高度
  treeheight_col = 15,  # 列聚类树高度（按样本表达聚类）
  fontsize = 12,  # 基础字体大小
  fontsize_row = 10,  # LR 对名称字体大小（避免重叠）
  fontsize_col = 12,  # 样本名字体大小
  color = heat_colors,  # 热图颜色梯度
  main = "LR Pairs Expression (DN2 vs NC)",  # 标题（明确对比组别）
  filename = "LR_pair_heatmap_DN2_vs_NC.pdf",  # 保存路径
  width = 12,  # 宽度（根据样本数量调整，避免过宽）
  height = ifelse(nrow(plot_data) > 20, 14, 10),  # 高度（根据 LR 对数量调整）
  border_color = NA,  # 隐藏单元格边框，更美观
  cluster_cols = TRUE,  # 样本按表达聚类
  cellwidth = 20,  # 单元格宽度
  cellheight = 12  # 单元格高度（避免 LR 对名称重叠）
)

# 13. 绘制热图（PNG格式，适合展示）
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
  dpi = 300  # 高清分辨率，适合展示
)

cat("热图已保存为 PDF 和 PNG 格式，仅包含 DN2 和 NC 两组样本！\n")
cat("保存文件：LR_pair_heatmap_DN2_vs_NC.pdf 和 LR_pair_heatmap_DN2_vs_NC.png\n")
