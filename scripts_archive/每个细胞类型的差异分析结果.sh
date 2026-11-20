# -------------------------- 加载必要的包 --------------------------
library(Seurat)
library(dplyr)
library(ClusterGVis)
library(org.Hs.eg.db)
library(ggsci)
library(ggplot2)
library(pheatmap)
library(tibble)

# -------------------------- 1. 通用配置 --------------------------
# 1.1 定义样本顺序和分组信息
sample_order <- c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6")

# 分组方案：DN1/DN2 vs NC
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

annotation_colors2 <- list(Group = c("DN1" = "#E41A1C", "DN2" = "#FF7F00", "NC" = "#377EB8"))
change_colors <- list(Change = c("up in case" = "#E41A1C", "down in case" = "#377EB8"))

# -------------------------- 2. 定义文件路径 --------------------------
file_paths <- c(
  "Activated_myofibroblasts" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Activated_myofibroblasts/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "B" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/B/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "Cortical_fibroblasts" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Cortical_fibroblasts/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "Distal_tubular_cell" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Distal_tubular_cell/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "M2_Macrophage" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/M2_Macrophage/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "Monocyte_cell" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Monocyte_cell/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "NK" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/NK/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "Plasma" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Plasma/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "Podocytes" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Podocytes/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "Proximal_tubular_cell" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/Proximal_tubular_cell/DE_list_Group2中DN比NC的差异_edgeR.csv",
  "T" = "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix/T/DE_list_Group2中DN比NC的差异_edgeR.csv"
)

# -------------------------- 3. 设置输出目录 --------------------------
out_path <- "/home/zhh2/YUBO-202508/Seurat_DEG_Analysis"
dir.create(out_path, showWarnings = FALSE)

# -------------------------- 4. 读取和预处理数据函数 --------------------------
read_and_preprocess_data <- function(file_path, cell_type) {
  cat("正在读取:", cell_type, "\n")
  
  # 读取CSV文件
  data <- read.csv(file_path, row.names = 1, check.names = FALSE)
  
  # 提取基因名
  genes <- rownames(data)
  
  # 提取表达量列（排除统计列）
  # 根据您的数据格式，列名格式如："D6_Proximal tubular cell"
  expr_cols <- colnames(data)[grepl("^[A-Z][0-9]+_", colnames(data))]
  
  if (length(expr_cols) == 0) {
    # 尝试另一种模式
    expr_cols <- colnames(data)[!grepl("logFC|pvalue|p.val|adj|padj|FDR", 
                                        colnames(data), ignore.case = TRUE)]
  }
  
  # 提取表达矩阵
  expr_matrix <- as.matrix(data[, expr_cols])
  rownames(expr_matrix) <- genes
  
  # 从列名中提取样本信息
  # 列名格式：D6_Proximal tubular cell
  sample_names <- sapply(colnames(expr_matrix), function(x) {
    # 提取下划线前的部分作为样本名
    strsplit(x, "_")[[1]][1]
  })
  
  # 简化列名
  colnames(expr_matrix) <- sample_names
  
  # 按照sample_order排序列
  expr_matrix <- expr_matrix[, intersect(sample_order, colnames(expr_matrix))]
  
  return(expr_matrix)
}

# -------------------------- 5. Seurat差异分析函数 --------------------------
perform_seurat_de_analysis <- function(expr_matrix, cell_type) {
  cat("\n=== 对", cell_type, "进行Seurat差异分析 ===\n")
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = expr_matrix,
    project = cell_type,
    min.cells = 1,
    min.features = 1
  )
  
  # 添加样本信息
  sample_names <- colnames(expr_matrix)
  seurat_obj$sample <- sample_names
  
  # 添加分组信息
  seurat_obj$group <- group_info2[sample_names, "Group"]
  
  # 合并DN1和DN2为DN组
  seurat_obj$dn_nc <- ifelse(seurat_obj$group %in% c("DN1", "DN2"), "DN", "NC")
  
  # 数据归一化
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # 数据缩放
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  
  # 设置分组标识
  Idents(seurat_obj) <- "dn_nc"
  
  # 使用Wilcoxon检验进行差异分析
  # 放宽阈值以检测更多差异基因
  de_markers <- FindMarkers(
    object = seurat_obj,
    ident.1 = "DN",
    ident.2 = "NC",
    test.use = "wilcox",
    logfc.threshold = 0.25,  # 降低logFC阈值
    min.pct = 0.1,          # 在至少10%的细胞中表达
    only.pos = FALSE,
    verbose = FALSE
  )
  
  # 添加基因名
  de_markers$gene <- rownames(de_markers)
  de_markers$cell_type <- cell_type
  
  # 使用更宽松的阈值标记上下调
  # 原阈值：log2FC > 1 且 p_val_adj < 0.05
  # 新阈值：log2FC > 0.5 且 p_val_adj < 0.1
  de_markers$change <- ifelse(
    de_markers$avg_log2FC > 0.5 & de_markers$p_val_adj < 0.1,
    "up in case",
    ifelse(
      de_markers$avg_log2FC < -0.5 & de_markers$p_val_adj < 0.1,
      "down in case",
      "not significant"
    )
  )
  
  # 提取显著差异基因（使用宽松阈值）
  sig_genes <- de_markers %>%
    filter(p_val_adj < 0.1 & abs(avg_log2FC) > 0.5)
  
  cat("  发现", nrow(sig_genes), "个显著差异基因 (log2FC|>0.5 & p.adj<0.1)\n")
  cat("  上调:", sum(sig_genes$avg_log2FC > 0), "个\n")
  cat("  下调:", sum(sig_genes$avg_log2FC < 0), "个\n")
  
  # 提取平均表达矩阵
  avg_exp <- AverageExpression(seurat_obj, assays = "RNA", slot = "data")[[1]]
  
  return(list(
    seurat_obj = seurat_obj,
    de_markers = de_markers,
    sig_genes = sig_genes,
    avg_expression = avg_exp
  ))
}

# -------------------------- 6. 基因聚类和GO富集函数 --------------------------
perform_clustering_go <- function(expr_matrix, cell_type, top_n = 20) {
  cat("\n=== 对", cell_type, "进行基因聚类和GO富集 ===\n")
  
  # 数据标准化（z-score）
  expr_scaled <- t(scale(t(expr_matrix)))
  
  # 检查数据维度
  if (nrow(expr_scaled) < 5) {
    cat("  警告：基因数量太少，跳过聚类分析\n")
    return(NULL)
  }
  
  # 设置聚类数量
  k <- min(4, nrow(expr_scaled) - 1)
  if (k < 2) k <- 2
  
  # 使用k-means聚类
  ck <- clusterData(
    obj = expr_scaled,
    clusterMethod = "kmeans",
    clusterNum = k
  )
  
  # 查看聚类结果
  cat("  聚类结果（各聚类基因数）：\n")
  print(table(ck$wide.res$cluster))
  
  # GO富集分析
  enrich_result <- tryCatch({
    enrichCluster(
      object = ck,
      OrgDb = org.Hs.eg.db,
      type = "BP",
      pvalueCutoff = 0.05,
      topn = 5,
      seed = 5201314
    )
  }, error = function(e) {
    cat("  GO富集分析错误:", e$message, "\n")
    return(NULL)
  })
  
  return(list(
    cluster_result = ck,
    go_enrichment = enrich_result
  ))
}

# -------------------------- 7. 可视化函数 --------------------------
visualize_heatmap <- function(expr_matrix, cell_type, change_info = NULL) {
  # 数据标准化
  expr_scaled <- t(scale(t(expr_matrix)))
  
  # 创建列注释（分组信息）
  sample_names <- colnames(expr_scaled)
  annotation_col <- data.frame(
    Group = group_info2[sample_names, "Group"],
    row.names = sample_names
  )
  
  # 如果有上下调信息，创建行注释
  annotation_row <- NULL
  if (!is.null(change_info)) {
    gene_changes <- change_info[rownames(expr_scaled), "change", drop = FALSE]
    colnames(gene_changes) <- "Change"
    annotation_row <- gene_changes
  }
  
  # 设置颜色
  annotation_colors <- list(
    Group = annotation_colors2$Group
  )
  
  if (!is.null(annotation_row)) {
    annotation_colors$Change <- change_colors$Change
  }
  
  # 绘制热图
  p <- pheatmap(
    expr_scaled,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 15,
    treeheight_col = 15,
    fontsize = 12,
    fontsize_row = 11,
    fontsize_col = 13,
    color = heat_colors,
    main = paste(cell_type, "Expression (DN vs NC)"),
    border_color = NA,
    cluster_cols = TRUE,
    cellwidth = 18,
    cellheight = 14,
    silent = FALSE
  )
  
  return(p)
}

visualize_cluster_go <- function(cluster_result, go_enrichment, cell_type) {
  # 准备termanno数据
  if (!is.null(go_enrichment)) {
    termanno <- data.frame(
      id = go_enrichment$group,
      term = go_enrichment$Description,
      pval = go_enrichment$pvalue,
      stringsAsFactors = FALSE
    )
  } else {
    termanno <- NULL
  }
  
  # 可视化1：热图和折线图
  pdf_file1 <- paste0(out_path, "/", cell_type, "_Cluster_Both.pdf")
  pdf(pdf_file1, height = 14, width = 18, onefile = FALSE)
  
  visCluster(
    object = cluster_result,
    plot.type = "both",
    column_names_rot = 45,
    markGenes.side = "left",
    line.side = "left",
    annoTerm.data = termanno,
    show_row_names = TRUE,
    go.col = rep(ggsci::pal_d3()(8), each = 4)[1:length(unique(go_enrichment$group))],
    go.size = 18
  )
  
  dev.off()
  cat("  聚类可视化1已保存:", pdf_file1, "\n")
  
  # 可视化2：只显示热图
  pdf_file2 <- paste0(out_path, "/", cell_type, "_Cluster_Heatmap.pdf")
  pdf(pdf_file2, height = 10, width = 12, onefile = FALSE)
  
  visCluster(
    object = cluster_result,
    plot.type = "heatmap",
    column_names_rot = 45,
    show_row_dend = FALSE,
    annoTerm.data = termanno,
    show_row_names = TRUE,
    line.side = "left",
    go.col = rep(ggsci::pal_d3()(8), each = 5),
    go.size = "pval"
  )
  
  dev.off()
  cat("  聚类可视化2已保存:", pdf_file2, "\n")
}

# -------------------------- 8. 主分析流程 --------------------------
analysis_results <- list()

# 先测试一个细胞类型
cell_type <- "Proximal_tubular_cell"  # 选择一个细胞类型测试

cat("\n", paste0(rep("=", 50), collapse = ""), "\n", sep = "")
cat("开始分析细胞类型:", cell_type, "\n")

# 1. 读取和预处理数据
expr_matrix <- read_and_preprocess_data(file_paths[cell_type], cell_type)

# 2. 保存原始表达矩阵
write.csv(expr_matrix, 
          file = paste0(out_path, "/", cell_type, "_expression_matrix.csv"))

# 3. Seurat差异分析
seurat_result <- perform_seurat_de_analysis(expr_matrix, cell_type)

if (!is.null(seurat_result)) {
  # 保存差异分析结果
  write.csv(seurat_result$de_markers,
            file = paste0(out_path, "/", cell_type, "_Seurat_DE_results.csv"))
  
  write.csv(seurat_result$sig_genes,
            file = paste0(out_path, "/", cell_type, "_Significant_DE_genes.csv"))
  
  # 4. 提取top差异基因进行聚类
  # 如果没有显著差异基因，使用所有基因的top 20（按p值排序）
  if (nrow(seurat_result$sig_genes) > 0) {
    top_genes <- seurat_result$sig_genes %>%
      arrange(p_val_adj) %>%
      head(20) %>%
      pull(gene)
  } else {
    # 如果没有显著差异基因，使用p值最小的20个基因
    cat("  没有显著差异基因，使用p值最小的20个基因\n")
    top_genes <- seurat_result$de_markers %>%
      arrange(p_val_adj) %>%
      head(20) %>%
      pull(gene)
  }
  
  if (length(top_genes) > 0) {
    top_expr <- expr_matrix[rownames(expr_matrix) %in% top_genes, , drop = FALSE]
    
    # 保存top基因表达矩阵
    write.csv(top_expr,
              file = paste0(out_path, "/", cell_type, "_Top20_genes_expression.csv"))
    
    # 5. 基因聚类和GO富集
    cluster_go_result <- perform_clustering_go(top_expr, cell_type)
    
    # 检查cluster_go_result是否存在
    if (!is.null(cluster_go_result)) {
      # 6. 可视化
      # 6.1 差异基因热图
      pdf_file <- paste0(out_path, "/", cell_type, "_DE_heatmap.pdf")
      pdf(pdf_file, width = 14, height = max(10, nrow(top_expr) * 0.5))
      
      # 设置行名为基因名
      change_info <- seurat_result$de_markers %>%
        filter(gene %in% top_genes) %>%
        select(gene, change)
      
      if (nrow(change_info) > 0) {
        rownames(change_info) <- change_info$gene
        change_info$gene <- NULL
      } else {
        change_info <- NULL
      }
      
      p <- visualize_heatmap(top_expr, cell_type, change_info)
      print(p)
      dev.off()
      cat("  差异基因热图已保存:", pdf_file, "\n")
      
      # 6.2 聚类和GO可视化
      if (!is.null(cluster_go_result$cluster_result)) {
        visualize_cluster_go(
          cluster_result = cluster_go_result$cluster_result,
          go_enrichment = cluster_go_result$go_enrichment,
          cell_type = cell_type
        )
        
        # 保存GO富集结果
        if (!is.null(cluster_go_result$go_enrichment)) {
          write.csv(cluster_go_result$go_enrichment,
                    file = paste0(out_path, "/", cell_type, "_GO_Enrichment.csv"))
        }
      }
    } else {
      cat("  聚类分析失败，跳过可视化\n")
    }
  }
  
  # 保存结果
  analysis_results[[cell_type]] <- list(
    expression_matrix = expr_matrix,
    seurat_analysis = seurat_result,
    clustering = if(exists("cluster_go_result")) cluster_go_result else NULL
  )
}

# -------------------------- 9. 继续分析其他细胞类型 --------------------------
# 分析其他细胞类型（跳过已分析的）
for (cell_type in names(file_paths)) {
  if (cell_type == "Proximal_tubular_cell") next  # 跳过已分析的
  
  cat("\n", paste0(rep("=", 50), collapse = ""), "\n", sep = "")
  cat("开始分析细胞类型:", cell_type, "\n")
  
  # 1. 读取和预处理数据
  expr_matrix <- read_and_preprocess_data(file_paths[cell_type], cell_type)
  
  # 2. 保存原始表达矩阵
  write.csv(expr_matrix, 
            file = paste0(out_path, "/", cell_type, "_expression_matrix.csv"))
  
  # 3. Seurat差异分析
  seurat_result <- perform_seurat_de_analysis(expr_matrix, cell_type)
  
  if (!is.null(seurat_result)) {
    # 保存差异分析结果
    write.csv(seurat_result$de_markers,
              file = paste0(out_path, "/", cell_type, "_Seurat_DE_results.csv"))
    
    write.csv(seurat_result$sig_genes,
              file = paste0(out_path, "/", cell_type, "_Significant_DE_genes.csv"))
    
    # 4. 提取top差异基因进行聚类
    # 如果没有显著差异基因，使用所有基因的top 20（按p值排序）
    if (nrow(seurat_result$sig_genes) > 0) {
      top_genes <- seurat_result$sig_genes %>%
        arrange(p_val_adj) %>%
        head(20) %>%
        pull(gene)
    } else {
      # 如果没有显著差异基因，使用p值最小的20个基因
      cat("  没有显著差异基因，使用p值最小的20个基因\n")
      top_genes <- seurat_result$de_markers %>%
        arrange(p_val_adj) %>%
        head(20) %>%
        pull(gene)
    }
    
    if (length(top_genes) > 0) {
      top_expr <- expr_matrix[rownames(expr_matrix) %in% top_genes, , drop = FALSE]
      
      # 保存top基因表达矩阵
      write.csv(top_expr,
                file = paste0(out_path, "/", cell_type, "_Top20_genes_expression.csv"))
      
      # 5. 基因聚类和GO富集
      cluster_go_result <- NULL
      tryCatch({
        cluster_go_result <- perform_clustering_go(top_expr, cell_type)
      }, error = function(e) {
        cat("  聚类分析失败:", e$message, "\n")
      })
      
      # 检查cluster_go_result是否存在
      if (!is.null(cluster_go_result)) {
        # 6. 可视化
        # 6.1 差异基因热图
        pdf_file <- paste0(out_path, "/", cell_type, "_DE_heatmap.pdf")
        pdf(pdf_file, width = 14, height = max(10, nrow(top_expr) * 0.5))
        
        # 设置行名为基因名
        change_info <- seurat_result$de_markers %>%
          filter(gene %in% top_genes) %>%
          select(gene, change)
        
        if (nrow(change_info) > 0) {
          rownames(change_info) <- change_info$gene
          change_info$gene <- NULL
        } else {
          change_info <- NULL
        }
        
        p <- visualize_heatmap(top_expr, cell_type, change_info)
        print(p)
        dev.off()
        cat("  差异基因热图已保存:", pdf_file, "\n")
        
        # 6.2 聚类和GO可视化
        if (!is.null(cluster_go_result$cluster_result)) {
          visualize_cluster_go(
            cluster_result = cluster_go_result$cluster_result,
            go_enrichment = cluster_go_result$go_enrichment,
            cell_type = cell_type
          )
          
          # 保存GO富集结果
          if (!is.null(cluster_go_result$go_enrichment)) {
            write.csv(cluster_go_result$go_enrichment,
                      file = paste0(out_path, "/", cell_type, "_GO_Enrichment.csv"))
          }
        }
      } else {
        cat("  聚类分析失败，跳过可视化\n")
      }
    }
    
    # 保存结果
    analysis_results[[cell_type]] <- list(
      expression_matrix = expr_matrix,
      seurat_analysis = seurat_result,
      clustering = cluster_go_result
    )
  }
}

# -------------------------- 10. 汇总分析结果 --------------------------
cat("\n", paste0(rep("=", 50), collapse = ""), "\n", sep = "")
cat("分析完成！生成汇总结果...\n")

# 创建汇总表
summary_list <- list()

for (cell_type in names(analysis_results)) {
  result <- analysis_results[[cell_type]]
  
  if (!is.null(result$seurat_analysis)) {
    sig_genes <- result$seurat_analysis$sig_genes
    de_markers <- result$seurat_analysis$de_markers
    
    # 使用不同阈值统计
    sig_strict <- de_markers %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1)
    sig_loose <- de_markers %>% filter(p_val_adj < 0.1 & abs(avg_log2FC) > 0.5)
    
    summary_list[[cell_type]] <- data.frame(
      Cell_Type = cell_type,
      Total_Genes = nrow(result$expression_matrix),
      Significant_Strict = nrow(sig_strict),  # 严格阈值
      Significant_Loose = nrow(sig_loose),    # 宽松阈值
      Upregulated = sum(sig_loose$avg_log2FC > 0),
      Downregulated = sum(sig_loose$avg_log2FC < 0),
      stringsAsFactors = FALSE
    )
  }
}

# 合并汇总表
if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  # 保存汇总结果
  write.csv(summary_df, 
            file = paste0(out_path, "/Analysis_Summary.csv"), 
            row.names = FALSE)
  
  cat("\n差异基因统计汇总:\n")
  print(summary_df)
  
  # 绘制汇总条形图
  if (require("reshape2", quietly = TRUE)) {
    library(reshape2)
    
    # 准备数据
    plot_data <- summary_df[, c("Cell_Type", "Upregulated", "Downregulated")]
    summary_melt <- melt(plot_data, id.vars = "Cell_Type")
    
    p_summary <- ggplot(summary_melt, aes(x = Cell_Type, y = value, fill = variable)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8")) +
      labs(title = "各细胞类型差异基因数量汇总 (DN vs NC)",
           subtitle = "阈值: |log2FC| > 0.5 & adj.p < 0.1",
           x = "细胞类型",
           y = "基因数量",
           fill = "变化方向") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(paste0(out_path, "/DE_Summary_Plot.pdf"), p_summary, width = 12, height = 8)
  }
}

# -------------------------- 11. 生成报告 --------------------------
cat("\n", paste0(rep("=", 50), collapse = ""), "\n", sep = "")
cat("分析报告\n")
cat("输出目录:", out_path, "\n")
cat("分析的细胞类型数量:", length(analysis_results), "\n")

cat("\n生成的文件:\n")
files <- list.files(out_path, pattern = "\\.(csv|pdf)$")
for (file in files) {
  cat("  -", file, "\n")
}

cat("\n分析完成！\n")
