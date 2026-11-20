# 1. 加载必备包
if (!require("pheatmap")) install.packages("pheatmap", dependencies = TRUE)
if (!require("RColorBrewer")) install.packages("RColorBrewer", dependencies = TRUE)
if (!require("stats")) install.packages("stats")
library(pheatmap)
library(RColorBrewer)
library(stats)
library(pheatmap)
library(RColorBrewer)
library(stats)
library(tibble) 

setwd("/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone")
# -------------------------- 2. 通用配置（所有文件共用） --------------------------
# 2.1 样本顺序和分组信息（与你提供的完全一致）
sample_order <- c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6")

# 分组方案1：DN vs NC（用于差异分析和注释）
group_info1 <- data.frame(
  Sample = sample_order,
  Group = c("DN", "DN", "DN", "DN", "DN", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = TRUE,
  row.names = sample_order
)

# 分组方案2：DN1/DN2 vs NC（仅用于注释）
group_info2 <- data.frame(
  Sample = sample_order,
  Group = c("DN2", "DN2", "DN1", "DN2", "DN1", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = TRUE,
  row.names = sample_order
)

# 2.2 颜色方案（与你提供的一致）

heat_colors <- colorRampPalette(c(
  "#006400", "#228B22", "#32CD32", "#90EE90",  # 下调：深绿→浅绿（对应down in case）
  "#F5F5F5",                                     # 中性色
  "#FFB6C1", "#FF6347", "#DC143C", "#B22222"    # 上调：浅红→深红（对应up in case）
))(100)

annotation_colors1 <- list(Group = c("DN" = "#E41A1C", "NC" = "#377EB8"))
annotation_colors2 <- list(Group = c("DN1" = "#E41A1C", "DN2" = "#FF7F00", "NC" = "#377EB8"))
change_colors <- list(Change = c("up in case" = "#E41A1C", "down in case" = "#377EB8"))

# 2.3 定义需要处理的3个文件路径
file_paths <- c(
  #"/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Signature/Matrix/1/Group2中DN比NC_Signatures.csv",
  "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Signature/Matrix/7/Group2中DN比NC_Signatures.csv",
  "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Signature/Matrix/9/Group2中DN比NC_Signatures.csv",
  "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Signature/Matrix/10/Group2中DN比NC_Signatures.csv"
)

# -------------------------- 3. 批量处理：差异筛选（前10上调+前10下调）+ 热图绘制 --------------------------
for (file in file_paths) {
  # 3.1 读取单个CSV文件并预处理
  # 提取文件编号（0/3/6）作为前缀，用于区分输出文件
  prefix <- paste0("Cell zone_", sub(".*Matrix/(\\d+)/.*", "\\1", file))
  # 读取数据（第一列作为行名，确保样本顺序与group_info匹配）
  data <- read.csv(
    file = file,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  
  # 数据检查：确保样本顺序与分组一致，缺失样本用NA填充
  missing_samples <- setdiff(sample_order, colnames(data))
  if (length(missing_samples) > 0) {
    warning(paste0("文件", file, "缺失样本：", paste(missing_samples, collapse = ","), "，已用NA填充"))
    na_df <- data.frame(matrix(NA, nrow = nrow(data), ncol = length(missing_samples), row.names = rownames(data)))
    colnames(na_df) <- missing_samples
    data <- cbind(data, na_df)
  }
  # 按sample_order重新排序列（确保与分组对应）
  data <- data[, sample_order, drop = FALSE]
  # 填充NA值（避免统计检验报错）
  if (any(is.na(data))) {
    data[is.na(data)] <- 0
    warning(paste0("文件", file, "存在NA值，已用0填充"))
  }
  
  # 3.2 显著差异筛选：先选P<0.05，再提前10上调+前10下调（合计20个）
  # 第一步：计算每个Signature的DN组vs NC组均值差异（判断上调/下调）和P值
  dn_cols <- which(group_info1$Group == "DN")
  nc_cols <- which(group_info1$Group == "NC")
  
  # 计算P值（t检验）
  p_values <- apply(data, 1, function(row) {
    t.test(row[dn_cols], row[nc_cols])$p.value
  })
  
  # 计算DN组vs NC组的均值差异（正=上调，负=下调）
  mean_diff <- apply(data, 1, function(row) {
    mean(row[dn_cols]) - mean(row[nc_cols])
  })
  
  # 第二步：整合数据，筛选P<0.05的显著Signature
  sig_full <- data.frame(
    Mean_Diff = mean_diff,
    P_value = p_values,
    row.names = rownames(data)
  ) %>%
    filter(P_value < 0.05) %>%  # 仅保留显著差异
    rownames_to_column("Signature")
  
  # 第三步：提取前10个上调和前10个下调（按均值差异绝对值排序，绝对值越大越显著）
  # 上调Signature（Mean_Diff正，取前10个）
  up_sig <- sig_full %>%
    filter(Mean_Diff > 0) %>%
    arrange(desc(abs(Mean_Diff))) %>%  # 按差异绝对值排序
    slice_head(n = 10)  # 取前10个最显著的
  
  # 下调Signature（Mean_Diff负，取前10个）
  down_sig <- sig_full %>%
    filter(Mean_Diff < 0) %>%
    arrange(desc(abs(Mean_Diff))) %>%
    slice_head(n = 10)
  
  # 合并为最终绘图数据（合计20个，不足则按实际数量）
  final_sig_names <- c(up_sig$Signature, down_sig$Signature)
  sig_data <- data[final_sig_names, , drop = FALSE]
  
  # 提示筛选结果
  cat(paste0("\n", prefix, " 筛选结果："), "\n")
  cat("  显著上调Signature（P<0.05）：", nrow(up_sig), "个（目标10个）\n")
  cat("  显著下调Signature（P<0.05）：", nrow(down_sig), "个（目标10个）\n")
  cat("  最终绘图Signature总数：", nrow(sig_data), "个\n")
  
  # 处理无显著Signature的极端情况
  if (nrow(sig_data) == 0) {
    warning(paste0("文件", file, "无P<0.05的显著Signature，使用全部数据绘图"))
    sig_data <- data
  }
  
  # 3.3 绘制热图（方案1：DN vs NC分组）
  # PNG格式
  png1 <- paste0(prefix, "_DNvsNC_Top20sig_heatmap.png")
  pheatmap(
    mat = sig_data,
    scale = "none",
    annotation_col = group_info1["Group", drop = FALSE],  # 列注释（样本分组）
    annotation_colors = annotation_colors1,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Top20 Signature (DN vs NC, P<0.05)"),  # 标题标注“Top20”
    filename = png1,
    width = 14,
    height = max(8, nrow(sig_data) * 0.6),  # 高度随Signature数量自适应（更紧凑）
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14,
    dpi = 300
  )
  
  # PDF格式
  pdf1 <- paste0(prefix, "_DNvsNC_Top20sig_heatmap.pdf")
  pheatmap(
    mat = sig_data,
    scale = "none",
    annotation_col = group_info1["Group", drop = FALSE],
    annotation_colors = annotation_colors1,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Top20 Signature (DN vs NC, P<0.05)"),
    filename = pdf1,
    width = 14,
    height = max(8, nrow(sig_data) * 0.6),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14
  )
  
  # 3.4 绘制热图（方案2：DN1/DN2 vs NC分组）
  # PNG格式
  png2 <- paste0(prefix, "_DN1DN2vsNC_Top20sig_heatmap.png")
  pheatmap(
    mat = sig_data,
    scale = "none",
    annotation_col = group_info2["Group", drop = FALSE],  # 列注释（细分DN1/DN2）
    annotation_colors = annotation_colors2,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Top20 Signature (DN1, DN2 vs NC, P<0.05)"),
    filename = png2,
    width = 14,
    height = max(8, nrow(sig_data) * 0.6),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14,
    dpi = 300
  )
  
  # PDF格式
  pdf2 <- paste0(prefix, "_DN1DN2vsNC_Top20sig_heatmap.pdf")
  pheatmap(
    mat = sig_data,
    scale = "none",
    annotation_col = group_info2["Group", drop = FALSE],
    annotation_colors = annotation_colors2,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(prefix, "Top20 Signature (DN1, DN2 vs NC, P<0.05)"),
    filename = pdf2,
    width = 14,
    height = max(8, nrow(sig_data) * 0.6),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14
  )
  
  # 输出当前文件的结果路径
  cat(paste0(prefix, " 结果路径：\n"))
  cat("  DN vs NC（Top20）：PNG=", normalizePath(png1), "，PDF=", normalizePath(pdf1), "\n")
  cat("  DN1/DN2 vs NC（Top20）：PNG=", normalizePath(png2), "，PDF=", normalizePath(pdf2), "\n\n")
}

cat("=== 所有文件处理完成！===\n")
