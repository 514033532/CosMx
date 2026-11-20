library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# Set working directory
setwd("/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_communication/stlearn/Group2_DN-VS-NC/lr_pair/")

# Read data
file_path <- "/home/zhh2/YUBO-202508/20250827/DN_share_results/20250826_results/groups/cell_communication/stlearn/Group2_DN-VS-NC/lr_pair/lr_value_DE_list_in_Group2_DN-VS-NC.csv"
data <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Define sample order
sample_order <- c("B1", "B3", "B4", "B5", "B6", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "A2", "D4", "D5", "D6")

# Create sample group data frames (English column names)
group_info1 <- data.frame(
  Sample = sample_order,
  Group = c("DN", "DN", "DN", "DN", "DN", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = FALSE,
  row.names = sample_order
)

group_info2 <- data.frame(
  Sample = sample_order,
  Group = c("DN2", "DN2", "DN1", "DN2", "DN1", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC", "NC"),
  stringsAsFactors = FALSE,
  row.names = sample_order
)

# Check required columns exist
required_cols <- c("Change", "Significance")
if (!all(required_cols %in% colnames(data))) {
  missing <- setdiff(required_cols, colnames(data))
  stop(paste("Required columns missing in data:", paste(missing, collapse = ", ")))
}

# Filter data with Significance = "Sig"
sig_data <- data %>% filter(Significance == "Sig")

# Filter 20 "up in case" and 20 "down in case" respectively
up_data <- sig_data %>% 
  filter(Change == "up in case") %>% 
  head(20)

down_data <- sig_data %>% 
  filter(Change == "down in case") %>% 
  head(20)

# Merge data
plot_data <- rbind(up_data, down_data)

# Check if there is enough data
if(nrow(plot_data) < 1) {
  stop("No qualified data available for heatmap plotting")
}

# Extract sample data and order by specified sample order
sample_data <- plot_data[, sample_order, drop = FALSE]

# Convert to matrix and standardize (row-wise)
matrix_data <- as.matrix(sample_data)
matrix_scaled <- t(scale(t(matrix_data)))  # Row-wise standardization

# Create row annotation (Change type of LR pairs)
annotation_row <- data.frame(
  Change = factor(plot_data$Change, levels = c("up in case", "down in case"))
)
rownames(annotation_row) <- rownames(plot_data)

# Set color schemes
heat_colors <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
annotation_colors1 <- list(
  Group = c("DN" = "#E41A1C", "NC" = "#377EB8")
)
annotation_colors2 <- list(
  Group = c("DN1" = "#E41A1C", "DN2" = "#FF7F00", "NC" = "#377EB8")
)

# Plot heatmap with group_info1 (DN vs NC) - Larger axis fonts
pheatmap(
  matrix_scaled,
  scale = "none",
  annotation_row = annotation_row,
  annotation_col = group_info1["Group", drop = FALSE],  # Use English column name
  annotation_colors = c(annotation_colors1, list(Change = c("up in case" = "#E41A1C", "down in case" = "#377EB8"))),
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 15,
  treeheight_col = 15,
  fontsize = 12,           # Global font size increased
  fontsize_row = 11,       # Y-axis (row names) font size - clearly visible
  fontsize_col = 13,       # X-axis (column names) font size - significantly increased
  color = heat_colors,
  main = "LR Pairs Expression (DN vs NC)",
  filename = "LR_pair_heatmap_group1.pdf",
  width = 14,              # Adjust width to fit large fonts
  height = 16,             # Adjust height to fit large fonts
  border_color = NA,
  cluster_cols = TRUE,     # Sample clustering
  cellwidth = 18,          # Increase cell width for better readability
  cellheight = 14          # Increase cell height for better readability
)

# Plot heatmap with group_info2 (DN1, DN2 vs NC) - Larger axis fonts
pheatmap(
  matrix_scaled,
  scale = "none",
  annotation_row = annotation_row,
  annotation_col = group_info2["Group", drop = FALSE],  # Use English column name
  annotation_colors = c(annotation_colors2, list(Change = c("up in case" = "#E41A1C", "down in case" = "#377EB8"))),
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 15,
  treeheight_col = 15,
  fontsize = 12,
  fontsize_row = 11,
  fontsize_col = 13,
  color = heat_colors,
  main = "LR Pairs Expression (DN1, DN2 vs NC)",
  filename = "LR_pair_heatmap_group2.pdf",
  width = 14,
  height = 16,
  border_color = NA,
  cluster_cols = TRUE,
  cellwidth = 18,
  cellheight = 14
)

# Save as PNG format (group_info1)
pheatmap(
  matrix_scaled,
  scale = "none",
  annotation_row = annotation_row,
  annotation_col = group_info1["Group", drop = FALSE],
  annotation_colors = c(annotation_colors1, list(Change = c("up in case" = "#E41A1C", "down in case" = "#377EB8"))),
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 15,
  treeheight_col = 15,
  fontsize = 12,
  fontsize_row = 11,
  fontsize_col = 13,
  color = heat_colors,
  main = "LR Pairs Expression (DN vs NC)",
  filename = "LR_pair_heatmap_group1.png",
  width = 14,
  height = 16,
  border_color = NA,
  cluster_cols = TRUE,
  cellwidth = 18,
  cellheight = 14,
  dpi = 300
)

# Save as PNG format (group_info2)
pheatmap(
  matrix_scaled,
  scale = "none",
  annotation_row = annotation_row,
  annotation_col = group_info2["Group", drop = FALSE],
  annotation_colors = c(annotation_colors2, list(Change = c("up in case" = "#E41A1C", "down in case" = "#377EB8"))),
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 15,
  treeheight_col = 15,
  fontsize = 12,
  fontsize_row = 11,
  fontsize_col = 13,
  color = heat_colors,
  main = "LR Pairs Expression (DN1, DN2 vs NC)",
  filename = "LR_pair_heatmap_group2.png",
  width = 14,
  height = 16,
  border_color = NA,
  cluster_cols = TRUE,
  cellwidth = 18,
  cellheight = 14,
  dpi = 300
)

cat("Heatmaps saved in PDF and PNG formats, including two grouping methods\n")
