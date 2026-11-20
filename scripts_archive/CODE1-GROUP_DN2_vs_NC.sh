# 加载所需库
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)


setwd("/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_type/proportion/Group1_DN2-VS-NC")
# 读取数据文件
file_path <- "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_type/proportion/Group1_DN2-VS-NC/CellType_prop_DE_list_in_Group1_DN2-VS-NC.csv"
data <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# 确保P值为数值型并处理可能的非数值值
data$Pvalue <- as.numeric(data$Pvalue)
if (any(is.na(data$Pvalue))) {
  warning("部分P值无法转换为数值型，已将其视为非显著")
  data$Pvalue[is.na(data$Pvalue)] <- 1  # 将NA值设为1，使其不被筛选
}

# 筛选出P值小于0.05的细胞类型
significant_cells <- rownames(data)[data$Pvalue < 0.05]
cat(paste0("发现", length(significant_cells), "个P值小于0.05的细胞类型\n"))

if (length(significant_cells) == 0) {
  cat("没有P值小于0.05的细胞类型，程序退出\n")
  quit()
}

# 创建样本分组数据框
group_info <- data.frame(
  样本 = c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6"),
  分组 = c("DN2", "DN2", "DN1", "DN2", "DN1", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = FALSE
)

# 定义需要分析的样本列
sample_columns <- c("A2", "B1", "B3", "B5", "D4", "D5", "D6", "F1", "F2", "F3", "F4", "F5", "F6", "F7")

# 检查样本列是否存在于数据中
missing_cols <- setdiff(sample_columns, colnames(data))
if (length(missing_cols) > 0) {
  warning(paste("数据中缺少以下样本列:", paste(missing_cols, collapse = ", ")))
  sample_columns <- intersect(sample_columns, colnames(data))
  if (length(sample_columns) == 0) {
    stop("没有找到有效的样本列，无法继续")
  }
}

# 提取样本数据并转换为长格式
sample_data <- data[, sample_columns, drop = FALSE] %>%
  rownames_to_column("CellType") %>%
  pivot_longer(cols = -CellType, names_to = "样本", values_to = "比例") %>%
  mutate(比例 = as.numeric(比例))  # 确保比例为数值型

# 合并分组信息并筛选需要比较的组别(仅保留NC和DN2)和显著的细胞类型
plot_data <- sample_data %>%
  merge(group_info, by = "样本") %>%
  filter(分组 %in% c("NC", "DN2") & CellType %in% significant_cells) %>%
  mutate(分组 = factor(分组, levels = c("NC", "DN2")))

# 检查是否有可用的绘图数据
if (nrow(plot_data) == 0) {
  stop("没有符合条件的绘图数据，请检查分组和样本匹配情况")
}

# 为每种显著的细胞类型创建独立的PDF和PNG文件
for (cell in significant_cells) {
  # 清理细胞类型名称中的特殊字符，用于文件名
  clean_cell_name <- gsub("[\\.| |/|\\\\]", "_", cell)
  
  # 获取当前细胞类型的数据和P值
  current_data <- plot_data %>% filter(CellType == cell)
  
  # 跳过没有数据的细胞类型
  if (nrow(current_data) == 0) {
    warning(paste("细胞类型", cell, "没有可用数据，已跳过"))
    next
  }
  
  p_value <- data[cell, "Pvalue"]
  p_label <- ifelse(p_value < 0.001, "P < 0.001", 
                   ifelse(p_value < 0.01, paste0("P = ", sprintf("%.3f", p_value)),
                          paste0("P = ", sprintf("%.4f", p_value))))
  
  # 创建图形
  p <- ggplot(current_data, aes(x = 分组, y = 比例, fill = 分组)) +
    # 箱线图
    geom_boxplot(width = 0.6, linewidth = 0.8, alpha = 0.7) +
    # 散点
    geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
    # 标题和坐标轴标签
    ggtitle(cell) +
    labs(x = "Group", y = "Cell proportion (%)") +
    # 填充颜色（NC和DN2的颜色定义）
    scale_fill_manual(values = c("NC" = "#377eb8", "DN2" = "#e41a1c")) +
    # 添加P值标注
    annotate("text", x = 1.5, y = max(current_data$比例) * 0.95, 
             label = p_label, size = 6, fontface = "bold") +
    # 主题设置（符合发表格式）
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, fontface = "bold"),
      axis.title.x = element_text(size = 20, margin = margin(t = 10), fontface = "bold"),
      axis.title.y = element_text(size = 20, margin = margin(r = 10), fontface = "bold"),
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    # 设置Y轴范围
    ylim(0, max(current_data$比例) * 1.1)
  
  # 保存为PDF和PNG格式
  pdf_file <- paste0(clean_cell_name, "_NC_vs_DN2.pdf")
  ggsave(pdf_file, plot = p, device = "pdf", 
         width = 12, height = 14, units = "cm", dpi = 300)
  
  png_file <- paste0(clean_cell_name, "_NC_vs_DN2.png")
  ggsave(png_file, plot = p, device = "png", 
         width = 12, height = 14, units = "cm", dpi = 300)
  
  cat(sprintf("已生成图形文件: %s 和 %s\n", pdf_file, png_file))
}

cat("所有显著差异的图形生成完成！\n")

