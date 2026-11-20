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
  "#006400", "#228B22", "#32CD32", "#90EE90",  # 下调：深绿→浅绿（对应down in case）
  "#F5F5F5",                                     # 中性色
  "#FFB6C1", "#FF6347", "#DC143C", "#B22222"    # 上调：浅红→深红（对应up in case）
))(100)
annotation_colors1 <- list(Group = c("DN" = "#E41A1C", "NC" = "#377EB8"))
annotation_colors2 <- list(Group = c("DN1" = "#E41A1C", "DN2" = "#FF7F00", "NC" = "#377EB8"))
change_colors <- list(Change = c("up in case" = "#E41A1C", "down in case" = "#377EB8"))


# 1.3 定义所有待分析的文件路径（共22个文件）
file_paths <- list(
  # gene_NMF相关
  gene_NMF_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/gene_NMF/Group2_DN-VS-NC/cNMF_DE_list_in_Group2_DN-VS-NC.csv",
  gene_NMF_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/gene_NMF/Group2_DN-VS-NC/cNMF_DE_heatmap_in_Group2_DN-VS-NC.csv",
  
  # spatial_distr/moran相关
  GeneMoran_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/spatial_distr/moran/Group2_DN-VS-NC/GeneMoran_DE_list_in_Group2_DN-VS-NC.csv",
  GeneMoran_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/spatial_distr/moran/Group2_DN-VS-NC/GeneMoran_DE_heatmap_in_Group2_DN-VS-NC.csv",
  
  # cell_colocate/celltype_enrich相关
  CellEnrich_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_colocate/celltype_enrich/Group2_DN-VS-NC/CellEnrich_DE_list_in_Group2_DN-VS-NC.csv",
  CellEnrich_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_colocate/celltype_enrich/Group2_DN-VS-NC/CellEnrich_DE_heatmap_in_Group2_DN-VS-NC.csv",
  
  # spatial_distr/centrality相关
  CellDegree_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/spatial_distr/centrality/Group2_DN-VS-NC/CellDegree_DE_list_in_Group2_DN-VS-NC.csv",
  CellDegree_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/spatial_distr/centrality/Group2_DN-VS-NC/CellDegree_DE_heatmap_in_Group2_DN-VS-NC.csv",
  CellCloseness_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/spatial_distr/centrality/Group2_DN-VS-NC/CellCloseness_DE_list_in_Group2_DN-VS-NC.csv",
  CellCloseness_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/spatial_distr/centrality/Group2_DN-VS-NC/CellCloseness_DE_heatmap_in_Group2_DN-VS-NC.csv",
  CellClustering_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/spatial_distr/centrality/Group2_DN-VS-NC/CellClustering_DE_list_in_Group2_DN-VS-NC.csv",
  CellClustering_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/spatial_distr/centrality/Group2_DN-VS-NC/CellClustering_DE_heatmap_in_Group2_DN-VS-NC.csv",
  
  # cell_colocate/zone_conn相关
  CellZone_conn_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_colocate/zone_conn/Group2_DN-VS-NC/CellZone_conn_DE_list_in_Group2_DN-VS-NC.csv",
  CellZone_conn_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_colocate/zone_conn/Group2_DN-VS-NC/CellZone_conn_DE_heatmap_in_Group2_DN-VS-NC.csv",
  
  # cell_type/proportion相关
  CellType_prop_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/proportion/Group2_DN-VS-NC/CellType_prop_DE_list_in_Group2_DN-VS-NC.csv",
  CellType_prop_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/proportion/Group2_DN-VS-NC/CellType_prop_DE_heatmap_in_Group2_DN-VS-NC.csv",
  
  # cell_zone/proportion相关
  CellZone_prop_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_zone/proportion/Group2_DN-VS-NC/CellZone_prop_DE_list_in_Group2_DN-VS-NC.csv",
  CellZone_prop_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_zone/proportion/Group2_DN-VS-NC/CellZone_prop_DE_heatmap_in_Group2_DN-VS-NC.csv",
  
  # cell_communication/cellchat相关
  lr_pair_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_communication/cellchat/Group2_DN-VS-NC/lr_pair/lr_value_DE_list_in_Group2_DN-VS-NC.csv",
  lr_pair_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_communication/cellchat/Group2_DN-VS-NC/lr_pair/lr_value_DE_heatmap_in_Group2_DN-VS-NC.csv",
  cell_pair_DE_list = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_communication/cellchat/Group2_DN-VS-NC/cell_pair/cell_lr_num_DE_list_in_Group2_DN-VS-NC.csv",
  cell_pair_DE_heatmap = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_communication/cellchat/Group2_DN-VS-NC/cell_pair/cell_lr_num_DE_heatmap_in_Group2_DN-VS-NC.csv"
)



# -------------------------- 2. 批量处理函数（每个文件调用一次） --------------------------
plot_heatmap_for_file <- function(file_path, prefix) {
  # 2.1 读取数据并检查格式
  cat(paste("开始处理文件：", file_path, "\n"))
  data <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  # 检查必需列（Change和Significance）
  required_cols <- c("Change", "Significance")
  if (!all(required_cols %in% colnames(data))) {
    missing <- setdiff(required_cols, colnames(data))
    stop(paste(prefix, "文件缺少必需列：", paste(missing, collapse = ", ")))
  }
  
  # 2.2 筛选显著差异数据（Significance = "Sig"）
  sig_data <- data %>% filter(Significance == "Sig")
  if (nrow(sig_data) == 0) {
    cat(paste(prefix, "文件无显著差异数据，跳过绘图\n"))
    return()
  }
  
  # 2.3 分别筛选前20个上调和下调数据（不足20个则取全部）
  up_data <- sig_data %>% filter(Change == "up in case") %>% head(20)
  down_data <- sig_data %>% filter(Change == "down in case") %>% head(20)
  plot_data <- rbind(up_data, down_data)
  
  if (nrow(plot_data) < 1) {
    cat(paste(prefix, "文件无符合条件的上调/下调数据，跳过绘图\n"))
    return()
  }
  
  # 2.4 数据预处理（提取样本列、转矩阵、行标准化）
  sample_data <- plot_data[, sample_order, drop = FALSE]
  matrix_data <- as.matrix(sample_data)
  matrix_scaled <- t(scale(t(matrix_data)))  # 行标准化（每行均值为0，标准差为1）
  
  # 2.5 定义行注释（上调/下调类型）
  annotation_row <- data.frame(
    Change = factor(plot_data$Change, levels = c("up in case", "down in case"))
  )
  rownames(annotation_row) <- rownames(plot_data)
  
  # 2.6 定义输出目录（与输入文件同目录）
  output_dir <- dirname(file_path)
  setwd(output_dir)
  
  # -------------------------- 3. 绘制热图（两种分组方案，PDF+PNG格式） --------------------------
  # 3.1 分组方案1：DN vs NC（PDF）
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_row = annotation_row,
    annotation_col = group_info1["Group", drop = FALSE],
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
    height = max(10, nrow(plot_data) * 0.5),  # 高度随数据行数自适应
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14
  )
  
  # 3.2 分组方案1：DN vs NC（PNG，300dpi）
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_row = annotation_row,
    annotation_col = group_info1["Group", drop = FALSE],
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
  
  # 3.3 分组方案2：DN1/DN2 vs NC（PDF）
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_row = annotation_row,
    annotation_col = group_info2["Group", drop = FALSE],
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
  
  # 3.4 分组方案2：DN1/DN2 vs NC（PNG，300dpi）
  pheatmap(
    matrix_scaled,
    scale = "none",
    annotation_row = annotation_row,
    annotation_col = group_info2["Group", drop = FALSE],
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

# -------------------------- 4. 批量执行所有文件 --------------------------
# 为每个文件定义唯一前缀（与分析类型对应）
prefixes <- names(file_paths)

# 循环调用函数处理每个文件
for (i in seq_along(file_paths)) {
  plot_heatmap_for_file(file_paths[[i]], prefixes[i])
}

cat("所有文件处理完成！\n")
