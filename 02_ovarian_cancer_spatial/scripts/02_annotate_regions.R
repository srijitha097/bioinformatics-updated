library(Seurat)
library(ggplot2)

ovary <- readRDS("results/ovarian_spatial_seurat.rds")

region_map <- c(
  "0" = "Macrophage / immune-rich",
  "1" = "Stromal / fibroblast-like",
  "2" = "Unresolved region",
  "3" = "Epithelial / tumor-like",
  "4" = "Epithelial / progenitor-like",
  "5" = "Stromal / fibroblast-rich"
)

# Convert cluster numbers to region labels
region_labels <- unname(
  region_map[as.character(ovary$seurat_clusters)]
)

# Add region annotation to each spatial spot
ovary$region <- region_labels

# Check annotations
print(table(ovary$region))

# Spatial annotation plot
p <- SpatialDimPlot(
  ovary,
  group.by = "region",
  label = TRUE,
  label.size = 3
) +
  ggtitle("Ovarian Cancer Spatial Regions")

ggsave(
  "results/06_spatial_annotated_regions.pdf",
  p,
  width = 10,
  height = 8
)

# Save annotated metadata
write.csv(
  ovary@meta.data,
  "results/spot_metadata_annotated.csv"
)

# Save annotated Seurat object
saveRDS(
  ovary,
  "results/ovarian_spatial_annotated.rds"
)

cat("\nAnnotation complete.\n")
