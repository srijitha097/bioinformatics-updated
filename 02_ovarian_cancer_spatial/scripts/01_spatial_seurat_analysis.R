suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(42)

matrix_dir <- "data/raw/filtered_feature_bc_matrix"
spatial_dir <- "data/raw/spatial"
results_dir <- "results"

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

required_matrix <- c(
  "barcodes.tsv.gz",
  "features.tsv.gz",
  "matrix.mtx.gz"
)

missing_matrix <- required_matrix[
  !file.exists(file.path(matrix_dir, required_matrix))
]

if (length(missing_matrix) > 0) {
  stop(
    "Missing matrix files: ",
    paste(missing_matrix, collapse = ", ")
  )
}

if (!file.exists(file.path(spatial_dir, "scalefactors_json.json"))) {
  stop("Missing spatial/scalefactors_json.json")
}

if (!file.exists(file.path(spatial_dir, "tissue_lowres_image.png"))) {
  stop("Missing spatial/tissue_lowres_image.png")
}

message("Reading Visium count matrix...")
counts <- Read10X(data.dir = matrix_dir)

if (is.list(counts)) {
  if ("Gene Expression" %in% names(counts)) {
    counts <- counts[["Gene Expression"]]
  } else {
    counts <- counts[[1]]
  }
}

ovary <- CreateSeuratObject(
  counts = counts,
  assay = "Spatial",
  project = "OvarianCancerSpatial",
  min.cells = 3,
  min.features = 200
)

message("Reading tissue image and spot coordinates...")
image <- Read10X_Image(
  image.dir = spatial_dir,
  image.name = "tissue_lowres_image.png",
  filter.matrix = TRUE
)

common_cells <- intersect(Cells(ovary), Cells(image))
ovary <- subset(ovary, cells = common_cells)
image <- image[common_cells]

ovary[["slice1"]] <- image

# Mitochondrial percentage
ovary[["percent.mt"]] <- PercentageFeatureSet(
  ovary,
  pattern = "^MT-"
)

# QC plots
pdf(
  file.path(results_dir, "01_qc_violin.pdf"),
  width = 10,
  height = 5
)
print(
  VlnPlot(
    ovary,
    features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  )
)
dev.off()

# Basic QC filtering.
# These are broad thresholds for removing very low-quality spots.
ovary <- subset(
  ovary,
  subset =
    nFeature_Spatial >= 400 &
    nCount_Spatial >= 500 &
    percent.mt < 25
)

message("Spots after QC: ", ncol(ovary))

# Standard Seurat preprocessing
ovary <- NormalizeData(
  ovary,
  assay = "Spatial",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

ovary <- FindVariableFeatures(
  ovary,
  assay = "Spatial",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

ovary <- ScaleData(
  ovary,
  assay = "Spatial",
  features = VariableFeatures(ovary),
  verbose = FALSE
)

ovary <- RunPCA(
  ovary,
  assay = "Spatial",
  features = VariableFeatures(ovary),
  npcs = 30,
  verbose = FALSE
)

ovary <- FindNeighbors(
  ovary,
  reduction = "pca",
  dims = 1:20,
  verbose = FALSE
)

ovary <- FindClusters(
  ovary,
  resolution = 0.5,
  verbose = FALSE
)

ovary <- RunUMAP(
  ovary,
  reduction = "pca",
  dims = 1:20,
  verbose = FALSE
)

# Spatial cluster plot
p_spatial <- SpatialDimPlot(
  ovary,
  images = "slice1",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 3
) +
  ggtitle("Ovarian cancer spatial clusters")

ggsave(
  file.path(results_dir, "02_spatial_clusters.pdf"),
  p_spatial,
  width = 9,
  height = 7
)

# UMAP
p_umap <- DimPlot(
  ovary,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Ovarian cancer Visium clusters")

ggsave(
  file.path(results_dir, "03_umap_clusters.pdf"),
  p_umap,
  width = 8,
  height = 6
)

# Cluster markers
DefaultAssay(ovary) <- "Spatial"

markers <- FindAllMarkers(
  ovary,
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25
)

write.csv(
  markers,
  file.path(results_dir, "cluster_markers.csv"),
  row.names = FALSE
)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  top10,
  file.path(results_dir, "top10_markers_per_cluster.csv"),
  row.names = FALSE
)

# Marker panel for broad tissue programs
marker_panel <- c(
  "EPCAM", "KRT8", "KRT18", "KRT19",
  "COL1A1", "COL1A2", "DCN", "COL3A1",
  "PTPRC", "CD3D", "CD68", "LST1",
  "PECAM1", "VWF", "KDR"
)

marker_panel <- marker_panel[
  marker_panel %in% rownames(ovary)
]

if (length(marker_panel) > 0) {
  p_markers <- SpatialFeaturePlot(
    ovary,
    features = marker_panel,
    images = "slice1",
    ncol = 3,
    alpha = c(0.2, 1)
  )

  ggsave(
    file.path(results_dir, "04_spatial_marker_panel.pdf"),
    p_markers,
    width = 12,
    height = ceiling(length(marker_panel) / 3) * 3.5
  )
}

# Spatially variable genes with Moran's I
spatial_features <- head(
  VariableFeatures(ovary),
  1000
)

ovary <- FindSpatiallyVariableFeatures(
  ovary,
  assay = "Spatial",
  features = spatial_features,
  selection.method = "moransi",
  verbose = FALSE
)

svg <- SpatiallyVariableFeatures(
  ovary,
  method = "moransi"
)

svg_table <- data.frame(
  gene = svg,
  stringsAsFactors = FALSE
)

write.csv(
  svg_table,
  file.path(results_dir, "spatially_variable_genes.csv"),
  row.names = FALSE
)

top_spatial <- head(svg, 6)

if (length(top_spatial) > 0) {
  p_svg <- SpatialFeaturePlot(
    ovary,
    features = top_spatial,
    images = "slice1",
    ncol = 2,
    alpha = c(0.2, 1)
  )

  ggsave(
    file.path(results_dir, "05_top_spatial_genes.pdf"),
    p_svg,
    width = 9,
    height = 10
  )
}

# Save metadata and object
write.csv(
  ovary@meta.data,
  file.path(results_dir, "spot_metadata.csv")
)

saveRDS(
  ovary,
  file.path(results_dir, "ovarian_spatial_seurat.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(results_dir, "sessionInfo.txt")
)

message("Analysis complete.")
message("Clusters: ", length(levels(Idents(ovary))))
message("Spots analyzed: ", ncol(ovary))
message("Review cluster markers before assigning biological region labels.")
