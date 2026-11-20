

# 强制指定函数归属，避免冲突
library(Seurat)
library(ClusterGVis)
library(org.Hs.eg.db)
library(ggsci)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
select <- dplyr::select  # 锁定dplyr的select

# ===================== 2. 通用配置 =====================
sample_order <- c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6")
group_info2 <- data.frame(
  Sample = sample_order,
  Group = c("DN2", "DN2", "DN1", "DN2", "DN1", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = FALSE,
  row.names = sample_order
)

heat_colors <- colorRampPalette(c(
  "#006400", "#228B22", "#32CD32", "#90EE90",
  "#F5F5F5",
  "#FFB6C1", "#FF6347", "#DC143C", "#B22222"
))(100)
annotation_colors2 <- list(Group = c("DN1" = "#E41A1C", "DN2" = "#FF7F00", "NC" = "#377EB8"))

# 待处理文件路径
file_paths <- c(
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Activated_myofibroblasts/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/B/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Cortical_fibroblasts/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Distal_tubular_cell/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/M2_Macrophage/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Monocyte_cell/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/NK/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Plasma/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Podocytes/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Proximal_tubular_cell/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/T/DE_list_Group2中DN比NC的差异_edgeR.csv"
)

# ===================== 3. 鲁棒的分析函数（全场景兼容） =====================
analyze_cell_type <- function(file_path) {
  tryCatch({
    # 基础信息提取
    cell_type <- stringr::str_split(file_path, "/")[[1]][length(stringr::str_split(file_path, "/")[[1]])-1]
    out_prefix <- paste0(cell_type, "_DN_vs_NC")
    out_dir <- dirname(file_path)
    cat("🔹 开始处理", cell_type, "\n")
    
    # Step 1: 读取并清洗数据（极致鲁棒）
    df <- read.csv(file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    if (nrow(df) == 0) {
      cat("❌ ", cell_type, "：文件为空，跳过\n\n")
      return(NULL)
    }
    # 统一基因列名
    gene_col <- ifelse("Gene" %in% colnames(df), "Gene", colnames(df)[1])
    colnames(df)[colnames(df) == gene_col] <- "Gene"
    # 提取样本名（兼容任意列名格式）
    sample_names <- stringr::str_extract(colnames(df)[-1], "^[^_]+")
    colnames(df)[-1] <- sample_names
    
    # Step 2: 筛选有效样本列（仅保留存在的列，避免all_of报错）
    valid_samples <- intersect(sample_names, sample_order)
    if (length(valid_samples) == 0) {
      cat("❌ ", cell_type, "：无匹配的样本列，跳过\n\n")
      return(NULL)
    }
    # 清洗数据：仅保留有效列+非空基因
    df_clean <- df %>%
      dplyr::select(Gene, dplyr::all_of(valid_samples)) %>%  # 仅选有效列
      dplyr::filter(!is.na(Gene) & Gene != "" & !duplicated(Gene)) %>%
      tidyr::replace_na(replace = setNames(rep(0, length(valid_samples)), valid_samples))  # NA填0
    if (nrow(df_clean) < 2) {
      cat("❌ ", cell_type, "：有效基因数＜2，跳过\n\n")
      return(NULL)
    }
    
    # Step 3: Seurat差异分析（放宽阈值，避免无差异基因）
    exp_matrix <- df_clean %>% tibble::column_to_rownames("Gene") %>% as.matrix()
    seurat_obj <- Seurat::CreateSeuratObject(
      counts = exp_matrix,
      meta.data = group_info2[valid_samples, , drop = FALSE],  # 仅关联有效样本
      project = cell_type,
      min.cells = 0, min.features = 0  # 放宽过滤
    )
    seurat_obj <- Seurat::NormalizeData(seurat_obj)
    Seurat::Idents(seurat_obj) <- "Group"
    
    # 差异分析（放宽阈值，确保有结果）
    deg_results <- Seurat::FindMarkers(
      object = seurat_obj,
      ident.1 = intersect(c("DN1", "DN2"), unique(seurat_obj$Group)),
      ident.2 = intersect("NC", unique(seurat_obj$Group)),
      test.use = "wilcox",
      logfc.threshold = 0.1,    # 降低FC阈值
      min.pct = 0.05,           # 降低表达量阈值
      only.pos = FALSE,
      pval.cutoff = 0.1         # 放宽p值阈值
    )
    # 保存差异分析结果（即使无显著基因）
    write.csv(deg_results, paste0(out_dir, "/", out_prefix, "_DEG_Seurat.csv"), row.names = TRUE)
    cat("✅ ", cell_type, "：差异分析结果已保存\n")
    
    # Step 4: 提取差异基因（兼容无显著基因的情况）
    if (nrow(deg_results) == 0) {
      top_deg_genes <- head(rownames(exp_matrix), 20)  # 无差异基因时取前20个表达基因
      cat("⚠️ ", cell_type, "：无显著差异基因，使用前20个表达基因\n")
    } else {
      # 按p值排序，取前50（不足则取全部）
      top_deg_genes <- deg_results %>%
        dplyr::arrange(p_val_adj) %>%
        head(min(50, nrow(deg_results))) %>%
        rownames()
    }
    if (length(top_deg_genes) < 2) {
      cat("❌ ", cell_type, "：可用基因数＜2，跳过聚类/富集\n\n")
      return(NULL)
    }
    
    # Step 5: 分组平均表达矩阵
    group_avg <- df_clean %>%
      dplyr::filter(Gene %in% top_deg_genes) %>%
      tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
      dplyr::left_join(group_info2[, c("Sample", "Group")], by = "Sample") %>%
      dplyr::filter(!is.na(Group)) %>%  # 过滤无分组的样本
      dplyr::group_by(Gene, Group) %>%
      dplyr::summarise(Mean_Expression = mean(Expression), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = "Group", values_from = "Mean_Expression") %>%
      tidyr::replace_na(replace = setNames(rep(0, length(unique(.$Group))), unique(.$Group)))
    write.csv(group_avg, paste0(out_dir, "/", out_prefix, "_Group_Avg_Expression.csv"), row.names = FALSE)
    cat("✅ ", cell_type, "：平均表达矩阵已保存\n")
    
    # Step 6: k-means聚类（前置检查，避免NULL）
    cluster_matrix <- group_avg %>% 
      tibble::column_to_rownames("Gene") %>% 
      as.matrix()
    if (nrow(cluster_matrix) < 2 | ncol(cluster_matrix) < 2) {
      cat("⚠️ ", cell_type, "：聚类矩阵维度不足，跳过聚类/富集\n\n")
      return(NULL)
    }
    # 执行聚类
    ck <- ClusterGVis::clusterData(
      obj = cluster_matrix,
      clusterMethod = "kmeans",
      clusterNum = min(4, nrow(cluster_matrix)),  # 聚类数≤基因数
      scale = TRUE
    )
    write.csv(ck$wide.res, paste0(out_dir, "/", out_prefix, "_Kmeans_Cluster.csv"), row.names = TRUE)
    cat("✅ ", cell_type, "：聚类结果已保存\n")
    
    # Step 7: GO富集（兼容基因注释失败）
    enrich_result <- tryCatch({
      ClusterGVis::enrichCluster(
        object = ck,
        OrgDb = org.Hs.eg.db,
        type = "BP",
        pvalueCutoff = 0.1,  # 放宽p值
        topn = 3,            # 减少topn
        seed = 5201314
      )
    }, error = function(e) {
      cat("⚠️ ", cell_type, "：GO富集失败：", e$message, "\n")
      return(NULL)
    })
    if (!is.null(enrich_result) & nrow(enrich_result) > 0) {
      termanno <- enrich_result %>%
        as.data.frame() %>%
        dplyr::select(group, Description, pvalue) %>%
        dplyr::rename(id = group, term = Description, pval = pvalue)
      write.csv(termanno, paste0(out_dir, "/", out_prefix, "_GO_BP_Enrichment.csv"), row.names = FALSE)
      cat("✅ ", cell_type, "：GO富集结果已保存\n")
    } else {
      termanno <- NULL
      cat("⚠️ ", cell_type, "：无有效GO富集结果\n")
    }
    
    # Step 8: 可视化（兼容无富集结果）
    vis_pdf <- paste0(out_dir, "/", out_prefix, "_Cluster_GO_Vis.pdf")
    pdf(vis_pdf, height = 12, width = 10)
    ClusterGVis::visCluster(
      object = ck,
      plot.type = "heatmap",  # 无富集时仅画热图
      column_names_rot = 45,
      show_row_dend = FALSE,
      markGenes.side = "left",
      annoTerm.data = if (!is.null(termanno)) termanno else NULL,
      line.side = "left",
      go.col = rep(ggsci::pal_d3()(8), each = 4)[1:if (!is.null(termanno)) length(unique(termanno$id)) else 1],
      go.size = "pval",
      heatColList = list(
        col_range = c(-2, 0, 2),
        col_color = c("#006400", "#F5F5F5", "#B22222")
      ),
      main = paste(cell_type, "DN vs NC Cluster")
    )
    dev.off()
    
    # Step 9: 差异基因热图
    plot_data <- df_clean %>%
      dplyr::filter(Gene %in% top_deg_genes) %>%
      tibble::column_to_rownames("Gene") %>%
      as.matrix()
    heatmap_file <- paste0(out_dir, "/", out_prefix, "_DEG_Heatmap.pdf")
    pheatmap::pheatmap(
      plot_data,
      annotation_col = group_info2[colnames(plot_data), , drop = FALSE],
      annotation_colors = annotation_colors2,
      show_rownames = TRUE,
      show_colnames = TRUE,
      treeheight_row = 15,
      treeheight_col = 15,
      fontsize = 12,
      fontsize_row = 11,
      fontsize_col = 13,
      color = heat_colors,
      main = paste(cell_type, "DEG Expression (DN vs NC)"),
      filename = heatmap_file,
      width = 14,
      height = max(8, nrow(plot_data) * 0.4),  # 降低高度系数，避免过高
      border_color = NA,
      cluster_cols = TRUE,
      cellwidth = 18,
      cellheight = 14
    )
    cat("✅ ", cell_type, "：可视化文件已保存\n\n")
    
  }, error = function(e) {
    cat("❌ 处理", cell_type, "失败：", substr(e$message, 1, 100), "\n\n")
  })
}

# ===================== 4. 批量执行 =====================
cat("=============== 开始批量分析 ===============\n")
for (file in file_paths) {
  if (file.exists(file)) {
    analyze_cell_type(file)
  } else {
    cat("❌ 文件不存在：", file, "\n\n")
  }
}
cat("=============== 所有文件处理完成 ===============\n")
cat("🎉 结果均保存至对应文件目录！")
