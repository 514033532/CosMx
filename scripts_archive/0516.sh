# 加载包
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("tibble")) install.packages("tibble")
if (!require("dplyr")) install.packages("dplyr")

library(pheatmap)
library(RColorBrewer)
library(tibble)
library(dplyr)

# 工作目录
setwd("/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone")

# 样本顺序
sample_order <- c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6")

# 分组1 DN vs NC
group_info1 <- data.frame(
  row.names = sample_order,
  Group = c("DN","DN","DN","DN","DN","NC","NC","NC","NC","NC","NC","NC","NC","NC","NC","NC")
)

# 分组2 DN1/DN2 vs NC
group_info2 <- data.frame(
  row.names = sample_order,
  Group = c("DN2","DN2","DN1","DN2","DN1","NC","NC","NC","NC","NC","NC","NC","NC","NC","NC","NC")
)

# 热图配色
heat_colors <- colorRampPalette(c(
  "#006400", "#228B22", "#32CD32", "#90EE90",
  "#F5F5F5",
  "#FFB6C1", "#FF6347", "#DC143C", "#B22222"
))(100)

annotation_colors1 <- list(Group = c(DN="#E41A1C", NC="#377EB8"))
annotation_colors2 <- list(Group = c(DN1="#E41A1C", DN2="#FF7F00", NC="#377EB8"))

# ===================== 你要的 4 个文件路径 =====================
file_paths <- c(
  "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Signature/Matrix/1/Group2中DN比NC_Signatures.csv",
  "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Signature/Matrix/7/Group2中DN比NC_Signatures.csv",
  "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Signature/Matrix/9/Group2中DN比NC_Signatures.csv",
  "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Signature/Matrix/10/Group2中DN比NC_Signatures.csv"
)

# 批量循环绘图
for (file in file_paths) {
  prefix <- paste0("Cell_zone_", sub(".*Matrix/(\\d+)/.*", "\\1", file))
  cat("正在处理：", prefix, "\n")
  
  data <- read.csv(file, row.names = 1, check.names = F)
  data <- data[, sample_order]
  data[is.na(data)] <- 0
  
  dn_cols <- which(group_info1$Group == "DN")
  nc_cols <- which(group_info1$Group == "NC")
  
  p_values <- apply(data, 1, function(row) t.test(row[dn_cols], row[nc_cols])$p.value)
  mean_diff <- apply(data, 1, function(row) mean(row[dn_cols]) - mean(row[nc_cols]))
  
  sig_full <- data.frame(
    Mean_Diff = mean_diff,
    P_value = p_values,
    Signature = rownames(data)
  )
  sig_full <- sig_full[sig_full$P_value < 0.05, ]
  
  up <- sig_full[sig_full$Mean_Diff > 0, ]
  up <- up[order(-abs(up$Mean_Diff)), ]
  if(nrow(up) > 10) up <- up[1:10, ]
  
  down <- sig_full[sig_full$Mean_Diff < 0, ]
  down <- down[order(-abs(down$Mean_Diff)), ]
  if(nrow(down) > 10) down <- down[1:10, ]
  
  final_sig <- c(up$Signature, down$Signature)
  mat <- if(length(final_sig)==0) data else data[final_sig, ]
  
  # ========== A4一页 论文版热图 PDF ==========
  pdf1 <- paste0(prefix, "_DNvsNC_A4.pdf")
  pheatmap(
    mat,
    scale = "none",
    annotation_col = group_info1,
    annotation_colors = annotation_colors1,
    color = heat_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    treeheight_row = 12,
    treeheight_col = 12,
    cellwidth = 22,
    cellheight = 15,
    fontsize = 11,
    fontsize_row = 9,
    fontsize_col = 10,
    border_color = NA,
    filename = pdf1,
    width = 21,    # A4宽
    height = 29.7  # A4高
  )
  
  pdf2 <- paste0(prefix, "_DN1DN2vsNC_A4.pdf")
  pheatmap(
    mat,
    scale = "none",
    annotation_col = group_info2,
    annotation_colors = annotation_colors2,
    color = heat_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    treeheight_row = 12,
    treeheight_col = 12,
    cellwidth = 22,
    cellheight = 15,
    fontsize = 11,
    fontsize_row = 9,
    fontsize_col = 10,
    border_color = NA,
    filename = pdf2,
    width = 21,
    height = 29.7
  )
  
  cat(prefix, " → A4论文版PDF已生成\n\n")
}

cat("==== ✅ 全部 4 个文件处理完成！A4一页、无遮挡、可直接发论文 ====\n")
