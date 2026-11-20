# 加载必备包
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)
if (!require("stringr")) install.packages("stringr", dependencies = TRUE)
if (!require("patchwork")) install.packages("patchwork", dependencies = TRUE)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)

# -------------------------- 1. 通用配置与增强版NA处理函数 --------------------------
file_paths <- c(
  "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Enrichment_results/0/GSEA_KEGG/Group2中DN比NC_GSEA_KEGG_enrichment.csv",
  "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Enrichment_results/3/GSEA_KEGG/Group2中DN比NC_GSEA_KEGG_enrichment.csv",
  "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_zone/Enrichment_results/6/GSEA_KEGG/Group2中DN比NC_GSEA_KEGG_enrichment.csv"
)

up_color <- "#E41A1C"
down_color <- "#377EB8"
fontsize <- 24

# 增强版NA处理函数
clean_na_data <- function(data, file_id) {
  data_char <- data %>% mutate(across(everything(), as.character))
  na_values <- c("", "NA", "N/A", "NULL", "na", "n/a", "null", "NaN", "Inf", "-Inf")
  na_summary <- data_char %>%
    summarise_all(~sum(. %in% na_values)) %>%
    tidyr::pivot_longer(everything(), names_to = "Column", values_to = "NA_Count") %>%
    filter(NA_Count > 0)
  
  if (nrow(na_summary) > 0) {
    warning(paste0("文件", file_id, "检测到NA值：\n", 
                  paste(apply(na_summary, 1, function(x) paste(x[1], ":", x[2], "个")), collapse = "\n")))
  }
  
  data_clean <- data %>%
    mutate(setSize = ifelse(is.na(setSize) | setSize %in% na_values, 
                           ifelse(all(is.na(setSize) | setSize %in% na_values), 
                                  20, 
                                  median(setSize[!is.na(setSize) & !setSize %in% na_values], na.rm = TRUE)), 
                           setSize)) %>%
    mutate(NES = ifelse(is.na(NES) | NES %in% na_values, 0, NES)) %>%
    mutate(p.adjust = ifelse(is.na(p.adjust) | p.adjust %in% na_values, 1, p.adjust)) %>%
    mutate(Description = ifelse(is.na(Description) | Description %in% na_values | Description == "", 
                               paste0("Unknown_Pathway_", row_number()), 
                               Description)) %>%
    filter(!(is.na(ID) | ID %in% na_values),
           !(is.na(enrichmentScore) | enrichmentScore %in% na_values),
           !(is.na(qvalue) | qvalue %in% na_values)) %>%
    mutate(across(c(setSize, NES, p.adjust, enrichmentScore, qvalue), as.numeric))
  
  cat(paste0("文件", file_id, "数据清洗完成：原始", nrow(data), "行 → 清洗后", nrow(data_clean), "行\n"))
  return(data_clean)
}

# -------------------------- 2. 批量处理每个文件（彻底修复bind_rows错误） --------------------------
# 通路名称自动分行函数（保持不变）
wrap_pathway_name <- function(name) {
  words <- str_split(name, " ")[[1]]
  if (length(words) <= 2) return(name)
  
  current_length <- 0
  wrapped_name <- ""
  
  for (word in words) {
    word_length <- nchar(word) + 1
    if (current_length + word_length > 25) {
      wrapped_name <- paste0(wrapped_name, "\n", word)
      current_length <- word_length - 1
    } else {
      wrapped_name <- paste0(wrapped_name, " ", word)
      current_length <- current_length + word_length
    }
  }
  
  str_trim(wrapped_name)
}

plot_list <- list()

for (i in seq_along(file_paths)) {
  file <- file_paths[i]
  file_id <- str_extract(file, "(?<=Enrichment_results/)[0-6](?=/GSEA_KEGG)")
  if (is.na(file_id)) file_id <- paste0("File_", i)
  
  # 读取数据并清洗NA
  tryCatch({
    enrich_raw <- read.csv(
      file = file,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA", "N/A", "NULL")
    )
    enrich_data <- clean_na_data(enrich_raw, file_id)
  }, error = function(e) {
    warning(paste0("文件", file_id, "读取失败：", e$message))
    enrich_data <- data.frame(
      ID = paste0("Error_", i, "_", 1:5),
      Description = paste0("Read_Error_Pathway_", 1:5),
      setSize = sample(10:30, 5),
      NES = sample(c(-2:-1, 1:2), 5, replace = TRUE),
      p.adjust = sample(0.01:0.05, 5, replace = TRUE)/100,
      enrichmentScore = sample(c(-2:-1, 1:2), 5, replace = TRUE),
      qvalue = sample(0.01:0.05, 5, replace = TRUE)/100,
      stringsAsFactors = FALSE
    )
  })
  
  # 筛选显著通路并分类
  enrich_processed <- enrich_data %>%
    mutate(Regulation = case_when(
      NES > 0 ~ "Up",
      NES < 0 ~ "Down",
      TRUE ~ "Neutral"
    )) %>%
    filter(p.adjust < 0.05, Regulation != "Neutral") %>%
    arrange(desc(abs(NES)))
  
  # -------------------------- 核心修复：处理上调/下调数据，确保字符型 --------------------------
  # 1. 处理上调数据（空数据框时手动指定列类型）
  if (nrow(filter(enrich_processed, Regulation == "Up")) == 0) {
    up_top5 <- data.frame(
      ID = character(),
      Description = character(),
      setSize = numeric(),
      NES = numeric(),
      p.adjust = numeric(),
      enrichmentScore = numeric(),
      qvalue = numeric(),
      Regulation = character(),
      Simplified_Desc = character(),  # 强制字符型
      stringsAsFactors = FALSE
    )
  } else {
    up_top5 <- enrich_processed %>%
      filter(Regulation == "Up") %>%
      slice_head(n = 5) %>%
      mutate(Simplified_Desc = str_remove(Description, "^hsa\\d+: ")) %>%
      mutate(Simplified_Desc = str_trunc(Simplified_Desc, 80, ellipsis = "...")) %>%
      # 用vapply强制返回字符型向量，彻底避免列表
      mutate(Simplified_Desc = vapply(Simplified_Desc, wrap_pathway_name, FUN.VALUE = character(1)))
  }
  
  # 2. 处理下调数据（同上，空数据框时强制字符型）
  if (nrow(filter(enrich_processed, Regulation == "Down")) == 0) {
    down_top5 <- data.frame(
      ID = character(),
      Description = character(),
      setSize = numeric(),
      NES = numeric(),
      p.adjust = numeric(),
      enrichmentScore = numeric(),
      qvalue = numeric(),
      Regulation = character(),
      Simplified_Desc = character(),  # 强制字符型
      stringsAsFactors = FALSE
    )
  } else {
    down_top5 <- enrich_processed %>%
      filter(Regulation == "Down") %>%
      slice_head(n = 5) %>%
      mutate(Simplified_Desc = str_remove(Description, "^hsa\\d+: ")) %>%
      mutate(Simplified_Desc = str_trunc(Simplified_Desc, 80, ellipsis = "...")) %>%
      mutate(Simplified_Desc = vapply(Simplified_Desc, wrap_pathway_name, FUN.VALUE = character(1)))
  }
  
  # 合并数据（此时两表Simplified_Desc均为字符型，无类型冲突）
  plot_data <- bind_rows(up_top5, down_top5) %>%
    mutate(
      Simplified_Desc = ifelse(duplicated(Simplified_Desc),
                              paste0(Simplified_Desc, " (", row_number(), ")"),
                              Simplified_Desc),
      Regulation = factor(Regulation, levels = unique(.$Regulation)),
      Simplified_Desc = factor(Simplified_Desc, levels = unique(.$Simplified_Desc))
    )
  
  # 处理无显著通路情况
  if (nrow(plot_data) == 0) {
    p <- ggplot() +
      geom_text(aes(x = 0, y = 0, label = "No significant\nUp/Down pathways"), 
                size = fontsize/3, color = "#666666", hjust = 0.5) +
      theme_void() +
      ggtitle(paste0("Cell zone ", file_id)) +
      theme(plot.title = element_text(size = fontsize, hjust = 0.5, margin = margin(b = 5)))
    plot_list[[i]] <- p
    next
  }
  
  # 绘制图形
  p <- ggplot(plot_data, aes(x = NES, y = Simplified_Desc)) +
    geom_point(aes(size = setSize, color = Regulation), shape = 16, stroke = 1) +
    scale_size_continuous(range = c(6, 12), name = "Gene Count") +
    scale_color_manual(values = c("Up" = up_color, "Down" = down_color), name = "Regulation") +
    scale_x_continuous(
      limits = c(min(-3, floor(min(plot_data$NES))), max(3, ceiling(max(plot_data$NES)))),
      breaks = seq(floor(min(plot_data$NES)), ceiling(max(plot_data$NES)), by = 1)
    ) +
    facet_wrap(~Regulation, scales = "free_y", ncol = 1, strip.position = "left") +
    labs(
      title = paste0("Cell zone ", file_id),
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_bw() +
    theme(
      plot.margin = margin(t = 10, b = 10, l = 60, r = 20, unit = "mm"),
      axis.text.y = element_text(
        size = fontsize , 
        hjust = 0, 
        lineheight = 1.1, 
        margin = margin(r = 10, unit = "mm")
      ),
      plot.title = element_text(size = fontsize, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      axis.title.x = element_text(size = fontsize , margin = margin(t = 8)),
      axis.text.x = element_text(size = fontsize ),
      strip.text.y = element_text(size = fontsize , face = "bold"),
      strip.background = element_rect(fill = "#F5F5F5"),
      legend.position = "bottom",
      legend.text = element_text(size = fontsize * 0.8),
      legend.title = element_text(size = fontsize * 0.9),
      legend.key.size = unit(1.2, "cm"),
      legend.box = "vertical",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.8, "lines")
    )
  
  plot_list[[i]] <- p
}

# -------------------------- 3. 组合图并保存 --------------------------
combined_plot <- wrap_plots(plot_list, nrow = 1, widths = rep(1, length(plot_list))) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "GSEA-KEGG Enrichment (Top 5 Up/Down Pathways, p.adjust < 0.05)",
    theme = theme(
      plot.title = element_text(size = fontsize * 1.1, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )

ggsave(
  filename = "KEGG_Enrichment_Wrapped_Name.png",
  plot = combined_plot,
  width = 35,
  height = 15,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "KEGG_Enrichment_Wrapped_Name.pdf",
  plot = combined_plot,
  width = 35,
  height = 15,
  bg = "white"
)

# -------------------------- 4. 输出结果提示 --------------------------
cat("\n=== 处理完成！===\n")
cat("1. 已彻底修复bind_rows()错误：\n")
cat("   - 用vapply强制返回字符型向量，替代易出问题的sapply\n")
cat("   - 空数据框时手动指定Simplified_Desc为字符型\n")
cat("2. 通路名称可自动分行，长名称完整显示\n")
cat("3. 生成文件：\n")
cat("   - PNG: ", normalizePath("KEGG_Enrichment_Wrapped_Name.png"), "\n")
cat("   - PDF: ", normalizePath("KEGG_Enrichment_Wrapped_Name.pdf"), "\n")
