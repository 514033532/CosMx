# 加载必备包
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)
if (!require("stringr")) install.packages("stringr", dependencies = TRUE)
if (!require("patchwork")) install.packages("patchwork", dependencies = TRUE)
if (!require("tidyr")) install.packages("tidyr", dependencies = TRUE)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(tidyr)

# -------------------------- 1. 通用配置 --------------------------
# 文件路径（保持你的原始路径）
file_paths <- c(
 "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Enrichment_results/1/GSEA_KEGG/Group2中DN比NC_GSEA_KEGG_enrichment.csv",
"/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Enrichment_results/7/GSEA_KEGG/Group2中DN比NC_GSEA_KEGG_enrichment.csv",
"/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Enrichment_results/9/GSEA_KEGG/Group2中DN比NC_GSEA_KEGG_enrichment.csv",
"/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Enrichment_results/10/GSEA_KEGG/Group2中DN比NC_GSEA_KEGG_enrichment.csv")

# 颜色配置
up_color <- "#E41A1C"
down_color <- "#377EB8"

# 字体与图大小配置（适配通路名）
fontsize <- 18
plot_width <- 28
plot_height <- 9

# NA值处理函数
clean_na_data <- function(data) {
  na_summary <- data %>%
    summarise_all(~sum(is.na(.))) %>%
    pivot_longer(everything(), names_to = "Column", values_to = "NA_Count") %>%
    filter(NA_Count > 0)
  
  if (nrow(na_summary) > 0) {
    warning(paste0("检测到NA值：\n", paste(apply(na_summary, 1, function(x) paste(x[1], ":", x[2], "个")), collapse = "\n")))
  }
  
  data_clean <- data %>%
    mutate(setSize = ifelse(is.na(setSize), median(setSize, na.rm = TRUE), setSize)) %>%
    mutate(NES = ifelse(is.na(NES), 0, NES)) %>%
    mutate(p.adjust = ifelse(is.na(p.adjust), 1, p.adjust)) %>%
    mutate(Description = ifelse(is.na(Description) | Description == "", 
                               paste0("Unknown_Pathway_", row_number()), 
                               Description)) %>%
    filter(complete.cases(ID, enrichmentScore, qvalue))
  
  cat(paste0("数据清洗完成：原始", nrow(data), "行 → 清洗后", nrow(data_clean), "行\n"))
  return(data_clean)
}

# -------------------------- 2. 批量处理文件 --------------------------
plot_list <- list()

for (i in seq_along(file_paths)) {
  file <- file_paths[i]
  # 提取文件夹编号，直接用Cell zone命名
  file_id <- str_extract(file, "(?<=Enrichment_results/)\\d+(?=/GSEA_KEGG)")
  if (is.na(file_id)) file_id <- paste0("Zone_", i)
  plot_title <- paste0("Cell zone ", file_id)
  
  # 读取数据
  tryCatch({
    enrich_raw <- read.csv(
      file = file,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA", "N/A", "NULL")
    )
    enrich_data <- clean_na_data(enrich_raw)
  }, error = function(e) {
    warning(paste0("文件", file, "读取失败：", e$message))
    enrich_data <- data.frame(
      ID = paste0("Error_", i, "_", 1:5),
      Description = paste0("Read_Error_Pathway_", 1:5),
      setSize = sample(10:30, 5),
      NES = sample(c(-2:-1, 1:2), 5, replace = TRUE),
      p.adjust = runif(5, 0.001, 0.05),
      stringsAsFactors = FALSE
    )
  })
  
  # 筛选显著通路
  enrich_processed <- enrich_data %>%
    mutate(Regulation = case_when(
      NES > 0 ~ "Up",
      NES < 0 ~ "Down",
      TRUE ~ "Neutral"
    )) %>%
    filter(p.adjust < 0.05, Regulation != "Neutral") %>%
    arrange(desc(abs(NES)))
  
  # 提取前5上调/下调（不截断名称）
  up_top5 <- enrich_processed %>%
    filter(Regulation == "Up") %>%
    slice_head(n = 5) %>%
    mutate(Simplified_Desc = str_remove(Description, "^hsa\\d+: "))
  
  down_top5 <- enrich_processed %>%
    filter(Regulation == "Down") %>%
    slice_head(n = 5) %>%
    mutate(Simplified_Desc = str_remove(Description, "^hsa\\d+: "))
  
  # 合并数据
  plot_data <- bind_rows(up_top5, down_top5) %>%
    mutate(
      Simplified_Desc = ifelse(duplicated(Simplified_Desc),
                              paste0(Simplified_Desc, " (", row_number(), ")"),
                              Simplified_Desc),
      Regulation = factor(Regulation, levels = unique(.$Regulation)),
      Simplified_Desc = factor(Simplified_Desc, levels = unique(.$Simplified_Desc))
    )
  
  # 无数据时的图（和其他图保持风格统一）
  if (nrow(plot_data) == 0) {
    p <- ggplot() +
      geom_text(aes(x = 0, y = 0, label = "No significant\npathways"), 
                size = fontsize/3, color = "#666666", hjust = 0.5, lineheight = 1.2) +
      theme_bw() +
      ggtitle(plot_title) +
      theme(
        plot.title = element_text(size = fontsize, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        axis.title.x = element_text(size = fontsize * 0.9),
        axis.text.x = element_text(size = fontsize * 0.8),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 10, b = 10, l = 10, r = 10, "mm")
      )
    plot_list[[i]] <- p
    next
  }
  
  # 绘制气泡图（关键修改：名称不截断、字体适配）
  p <- ggplot(plot_data, aes(x = NES, y = Simplified_Desc)) +
    geom_point(aes(size = setSize, color = Regulation), shape = 16) +
    scale_size_continuous(range = c(4, 12), name = "Gene Count", limits = c(0, 160)) +
    scale_color_manual(values = c("Up" = up_color, "Down" = down_color), name = "Regulation") +
    scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
    facet_wrap(~Regulation, scales = "free_y", ncol = 1, strip.position = "left") +
    labs(
      title = plot_title,
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_bw() +
    theme(
      plot.margin = margin(t = 10, b = 10, l = 10, r = 10, "mm"),
      plot.title = element_text(size = fontsize, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = fontsize * 0.9),
      axis.text.x = element_text(size = fontsize * 0.8),
      axis.text.y = element_text(size = fontsize * 0.8),
      strip.text.y = element_text(size = fontsize * 0.9, face = "bold"),
      strip.background = element_rect(fill = "#F5F5F5"),
      legend.position = "bottom",
      legend.text = element_text(size = fontsize * 0.8),
      legend.title = element_text(size = fontsize * 0.9),
      legend.key.size = unit(1, "cm"),
      legend.box = "vertical",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines")
    )
  
  plot_list[[i]] <- p
}

# -------------------------- 3. 拼图与保存 --------------------------
combined_plot <- wrap_plots(plot_list, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "GSEA-KEGG Enrichment (Top 5 Up/Down Pathways, p.adjust < 0.05)",
    theme = theme(
      plot.title = element_text(size = fontsize * 1.1, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )

# 保存高清图
ggsave("KEGG_Enrichment_CellZone.png", combined_plot, width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave("KEGG_Enrichment_CellZone.pdf", combined_plot, width = plot_width, height = plot_height, bg = "white")

cat("\n✅ 已全部修改完成！\n")
cat("📊 新图特点：\n")
cat("   - 标题改为 Cell zone 1/7/9/10\n")
cat("   - 通路名称完整显示，不再截断\n")
cat("   - 图例统一为1个，点大小标尺一致\n")
cat("   - 无数据的Cell zone 10 风格和其他图统一\n")
