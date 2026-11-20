# 加载所需库
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# 读取数据文件
file_path <- "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_type/proportion/Group1_DN1-VS-NC/CellType_prop_DE_list_in_Group1_DN1-VS-NC.csv"
data <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# 筛选出P值小于0.05的细胞类型
data$Pvalue <- as.numeric(data$Pvalue)  # 确保P值为数值型
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

# 提取样本数据并转换为长格式
sample_columns <- c("A2", "B4", "B6", "D4", "D5", "D6", "F1", "F2", "F3", "F4", "F5", "F6", "F7")
sample_data <- data[, sample_columns, drop = FALSE] %>%
  rownames_to_column("CellType") %>%
  pivot_longer(cols = -CellType, names_to = "样本", values_to = "比例")

# 合并分组信息并筛选需要比较的组别和显著的细胞类型
plot_data <- sample_data %>%
  merge(group_info, by = "样本") %>%
  filter(分组 %in% c("NC", "DN1") & CellType %in% significant_cells) %>%
  mutate(分组 = factor(分组, levels = c("NC", "DN1")))

# 为每种显著的细胞类型创建独立的PDF和PNG文件
for (cell in significant_cells) {
  # 清理细胞类型名称中的特殊字符，用于文件名
  clean_cell_name <- gsub("\\.| ", "_", cell)
  
  # 获取当前细胞类型的数据和P值
  current_data <- plot_data %>% filter(CellType == cell)
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
    # 填充颜色
    scale_fill_manual(values = c("NC" = "#377eb8", "DN1" = "#e41a1c")) +
    # 添加P值标注
    annotate("text", x = 1.5, y = max(current_data$比例) * 0.95, 
             label = p_label, size = 6, fontface = "bold") +
    # 主题设置（符合发表格式）
    theme(
      # 标题设置：20号
      plot.title = element_text(size = 20, hjust = 0.5, fontface = "bold"),
      # 坐标轴标题：20号
      axis.title.x = element_text(size = 20, margin = margin(t = 10), fontface = "bold"),
      axis.title.y = element_text(size = 20, margin = margin(r = 10), fontface = "bold"),
      # 坐标轴刻度标签：14号
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      # 坐标轴线条
      axis.line = element_line(color = "black", linewidth = 0.8),
      # 坐标轴刻度
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      # 封闭框架
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      # 隐藏图例
      legend.position = "none",
      # 去除网格线
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # 背景色
      panel.background = element_blank(),
      # 整体图形边距
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    # 设置Y轴范围
    ylim(0, max(current_data$比例) * 1.1)
  
  # 保存为PDF（适合印刷出版）
  pdf_file <- paste0(clean_cell_name, "_NC_vs_DN1.pdf")
  ggsave(pdf_file, plot = p, device = "pdf", 
         width = 12, height = 14, units = "cm", dpi = 300)
  
  # 保存为PNG（适合电子展示）
  png_file <- paste0(clean_cell_name, "_NC_vs_DN1.png")
  ggsave(png_file, plot = p, device = "png", 
         width = 12, height = 14, units = "cm", dpi = 300)
  
  cat(sprintf("已生成图形文件: %s 和 %s\n", pdf_file, png_file))
}

cat("所有显著差异的图形生成完成！\n")
    
