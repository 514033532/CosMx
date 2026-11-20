# 加载必需的包
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggrepel)

base_dir <- "/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Matrix"

cell_config <- data.frame(
  cell_name = c("B_cell", "NK_cell", "Mast_cell", "Proximal_tubular_cell",
                "Vascular_Endothelial_Cell", "Cortical_fibroblasts", "M2_Macrophage", "Monocyte_cell"),
  folder = c("B", "NK", "Mast", "Proximal_tubular_cell", "Vascular_Endothelial_Cell",
             "Cortical_fibroblasts", "M2_Macrophage", "Monocyte_cell"),
  stringsAsFactors = FALSE
)

get_colname <- function(df, patterns, exact = FALSE) {
  all_names <- names(df)
  for (pat in patterns) {
    idx <- if(exact) which(all_names == pat) else grep(pat, all_names, ignore.case = TRUE)
    if (length(idx) > 0) return(all_names[idx[1]])
  }
  return(NULL)
}

plot_volcano_labeled <- function(cell_name, folder, top_n_label = 5) {
  file_path <- file.path(base_dir, folder, "DE_list_Group2中DN比NC的差异_edgeR.csv")
  if (!file.exists(file_path)) {
    return(ggplot() + annotate("text", x=0, y=0, label=paste(cell_name, "\n文件不存在"), size=6) + theme_void())
  }
  
  deg <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  names(deg) <- tolower(trimws(names(deg)))
  
  gene_col <- get_colname(deg, c("genesymbol", "gene", "symbol"))
  if (is.null(gene_col)) deg$Gene <- rownames(deg) else deg$Gene <- deg[[gene_col]]
  
  logfc_col <- get_colname(deg, c("logfc", "log2fc", "logfoldchange"))
  if (is.null(logfc_col)) {
    return(ggplot() + annotate("text", x=0, y=0, label=paste(cell_name, "\n无 logFC 列"), size=6) + theme_void())
  }
  deg$logFC <- as.numeric(deg[[logfc_col]])
  
  fdr_col <- get_colname(deg, c("fdr", "adj.p.val", "padj", "bh"))
  if (is.null(fdr_col)) {
    pval_col <- get_colname(deg, c("pvalue", "p.value"))
    if (!is.null(pval_col)) {
      deg$FDR <- as.numeric(deg[[pval_col]])
    } else {
      return(ggplot() + annotate("text", x=0, y=0, label=paste(cell_name, "\n无FDR/P值"), size=6) + theme_void())
    }
  } else {
    deg$FDR <- as.numeric(deg[[fdr_col]])
  }
  
  deg <- deg[!is.na(deg$logFC) & !is.na(deg$FDR) & !is.na(deg$Gene), ]
  if(nrow(deg) == 0) {
    return(ggplot() + annotate("text", x=0, y=0, label=paste(cell_name, "\n无有效数据"), size=6) + theme_void())
  }
  
  deg$sig_group <- case_when(
    deg$FDR < 0.05 & deg$logFC > 1   ~ "Up",
    deg$FDR < 0.05 & deg$logFC < -1  ~ "Down",
    TRUE ~ "Not_Sig"
  )
  deg$log10FDR <- -log10(deg$FDR)
  deg$log10FDR[is.infinite(deg$log10FDR)] <- max(deg$log10FDR[is.finite(deg$log10FDR)], na.rm = TRUE)
  
  up_sig <- deg %>% filter(sig_group == "Up") %>% arrange(FDR) %>% head(top_n_label)
  down_sig <- deg %>% filter(sig_group == "Down") %>% arrange(FDR) %>% head(top_n_label)
  label_df <- bind_rows(up_sig, down_sig)
  
  # ========== 调整后的图形参数 ==========
  p <- ggplot(deg, aes(x = logFC, y = log10FDR, color = sig_group)) +
    geom_point(size = 1.5, alpha = 0.8) +  # 点变大、更不透明（更黑/更浓）
    scale_color_manual(values = c("Up" = "#D62828",        # 更深的红
                                  "Down" = "#023E8A",      # 更深的蓝
                                  "Not_Sig" = "#4A4A4A")) + # 深灰色（比之前更黑）
    geom_vline(xintercept = c(-1, 1), linetype = 2, color = "gray30", size = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "gray30", size = 0.6) +
    geom_text_repel(data = label_df, 
                    aes(x = logFC, y = log10FDR, label = Gene),
                    size = 3.5,          # 基因标签字体调大
                    box.padding = 0.4,
                    point.padding = 0.3,
                    max.overlaps = 12,
                    segment.color = "gray40",
                    segment.size = 0.4) +
    labs(title = gsub("_", " ", cell_name), 
         x = expression(log[2]("Fold Change")), 
         y = expression(-log[10]("FDR"))) +
    theme_bw(base_size = 12) +  # 基础字体从原来的8左右提高到12
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "none",
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 10),
          panel.grid.minor = element_blank())
  return(p)
}

# 批量生成火山图对象
plot_list <- list()
for (i in 1:nrow(cell_config)) {
  plot_list[[i]] <- plot_volcano_labeled(cell_config$cell_name[i], cell_config$folder[i], top_n_label = 5)
  cat(paste0("✅ ", cell_config$cell_name[i], "\n"))
}

# 公共图例（同样增大字体）
temp_plot <- ggplot(data.frame(x=0,y=0,group=factor(c("Up","Down","Not_Sig"))), 
                    aes(x,y,color=group)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("Up"="#D62828","Down"="#023E8A","Not_Sig"="#4A4A4A")) +
  theme_bw(base_size = 12) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11))
legend <- get_legend(temp_plot)

plot_list_no_leg <- lapply(plot_list, function(p) p + theme(legend.position = "none"))

combined <- wrap_plots(plot_list_no_leg, ncol = 4, nrow = 2) +
  plot_annotation(title = "DKD vs NC: Volcano plots across cell types (Top 5 genes labeled)",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))

final_plot <- plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.08))

# 保存（增大画布尺寸以容纳更大的字体）
ggsave("Combined_Volcano_Plots_large_font.pdf", final_plot, width = 16, height = 10, bg = "white")
ggsave("Combined_Volcano_Plots_large_font.png", final_plot, width = 16, height = 10, dpi = 300, bg = "white")

cat("\n🎉 已完成！字体调大、散点变大变深。输出文件：Combined_Volcano_Plots_large_font.pdf/png\n")
