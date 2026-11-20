

library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
setwd("/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_communication/stlearn/Group2_DN-VS-NC/lr_pair/")
# 读取数据
file_path <- "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_communication/stlearn/Group2_DN-VS-NC/lr_pair/lr_value_DE_list_in_Group2_DN-VS-NC.csv"

data <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# 确保P值为数值型并处理可能的非数值值
data$Pvalue <- as.numeric(data$Pvalue)
if (any(is.na(data$Pvalue))) {
  warning("部分P值无法转换为数值型，已将其视为非显著")
  data$Pvalue[is.na(data$Pvalue)] <- 1  # 将NA值设为1，使其不被筛选
}

# 筛选出P值小于0.05的细胞类型
significant_cells <- rownames(data)[data$Pvalue < 0.0005]
cat(paste0("发现", length(significant_cells), "个P值小于0.05的细胞类型\n"))

if (length(significant_cells) == 0) {
  cat("没有P值小于0.05的细胞类型，程序退出\n")
  quit()
}









significant_cells <-rownames(data)[1:15]
















# 创建样本分组数据框
group_info <- data.frame(
  样本 = c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6"),
  分组 = c("DN", "DN", "DN", "DN", "DN", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = FALSE
)

# 定义需要分析的样本列
sample_columns <- c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6")

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
  filter(分组 %in% c("NC", "DN") & CellType %in% significant_cells) %>%
  mutate(分组 = factor(分组, levels = c("NC", "DN")))

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
# 定义形状方案：NC组用圆点，DN组用三角形
shape_types <- c(
  "NC" = 16,  # 实心圆点
  "DN" = 17   # 实心三角形
)

# 定义颜色方案：深绿色、暗红色
group_colors <- c(
  "NC" = "#006400",    # 深绿色 (NC组)
  "DN" = "#8B0000"     # 暗红色 (DN组)
)

# 创建图形
p <- ggplot(current_data, aes(x = 分组, y = 比例, fill = 分组, shape = 分组)) +
  # 箱线图
  geom_boxplot(width = 0.6, linewidth = 0.8, alpha = 0.7) +
  # 散点 - 应用形状和颜色
  geom_jitter(aes(color = 分组), width = 0.1, size = 3, alpha = 0.8) +
  # 标题：超过8字符自动换行，按空格拆分（不截断单词）
  ggtitle(str_wrap(cell, width = 28)) +  # 关键：width=8，字符超8则自动换行
  # 坐标轴标签
  labs(x = "Group", y = "cell communication Score") +
  # 填充颜色（箱线图）
  scale_fill_manual(values = group_colors) +
  # 点的颜色
  scale_color_manual(values = group_colors) +
  # 点的形状
  scale_shape_manual(values = shape_types) +
  # 添加P值标注
  annotate("text", x = 1.5, y = max(current_data$比例) * 0.95, 
           label = p_label, size = 6, fontface = "bold") +
  # 主题设置：适配多行标题，避免拥挤
  theme(
    # 标题：增加行高和底部边距，适配换行后的显示
    plot.title = element_text(
      size = 20, 
      hjust = 0.5, 
      fontface = "bold",
      lineheight = 1.4,  # 增大行间距，防止多行标题上下重叠
      margin = margin(b = 20)  # 增加底部空白，避免标题压箱线图
    ),
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
    plot.margin = unit(c(1.5, 1, 1, 1), "cm")  # 顶部边距加宽，容纳多行标题
  ) +
  ylim(0, max(current_data$比例) * 1.1)
  
  # 保存为PDF和PNG格式
  pdf_file <- paste0(clean_cell_name, "_NC_vs_DN.pdf")
  ggsave(pdf_file, plot = p, device = "pdf", 
         width = 12, height = 14, units = "cm", dpi = 300)
  
  png_file <- paste0(clean_cell_name, "_NC_vs_DN.png")
  ggsave(png_file, plot = p, device = "png", 
         width = 12, height = 14, units = "cm", dpi = 300)
  
  cat(sprintf("已生成图形文件: %s 和 %s\n", pdf_file, png_file))
}

cat("所有显著差异的图形生成完成！\n")
    
