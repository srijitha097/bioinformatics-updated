suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(42)

input_dir <- "data/raw/filtered_feature_bc_matrix"
out_dir <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_dir)) {
  stop(
    "Missing 10x matrix directory: ", input_dir
  )
}

# -----------------------------
# 1. Load 10x count matrix
# -----------------------------
data_10x <- Read10X("data/raw/filtered_feature_bc_matrix")

if (is.list(data_10x)) {
  counts <- data_10x[["Gene Expression"]]
} else {
  counts <- data_10x
}
brain <- CreateSeuratObject(
  counts = counts,
  project = "E18_Mouse_Brain",
  min.cells = 3,
  min.features = 200
)
# Mouse mitochondrial genes begin with "mt-".
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-")

# -----------------------------
# 2. QC visualization
# -----------------------------
pdf(file.path(out_dir, "01_qc_violin.pdf"), width = 11, height = 4)
print(VlnPlot(
  brain,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.05
))
dev.off()

pdf(file.path(out_dir, "02_qc_scatter.pdf"), width = 10, height = 4)
p1 <- FeatureScatter(brain, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(brain, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p1 + p2)
dev.off()

# -----------------------------
# 3. Filter cells 
# Inspect QC distributions and tune for the selected dataset.
# -----------------------------
brain <- subset(
  brain,
  subset =
    nFeature_RNA >= 300 &
    nFeature_RNA <= 6500 &
    nCount_RNA >= 500 &
    percent.mt <= 15
)

# -----------------------------
# 4. Normalize + HVGs
# -----------------------------
brain <- NormalizeData(
  brain,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

brain <- FindVariableFeatures(
  brain,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# -----------------------------
# 5. Scale + PCA
# -----------------------------
brain <- ScaleData(brain, features = VariableFeatures(brain), verbose = FALSE)
brain <- RunPCA(brain, features = VariableFeatures(brain), npcs = 30, verbose = FALSE)

pdf(file.path(out_dir, "03_pca.pdf"), width = 6, height = 5)
print(DimPlot(brain, reduction = "pca"))
dev.off()

pdf(file.path(out_dir, "04_elbow.pdf"), width = 6, height = 5)
print(ElbowPlot(brain, ndims = 30))
dev.off()

# -----------------------------
# 6. Graph + clustering + UMAP
# -----------------------------
dims_use <- 1:20
brain <- FindNeighbors(brain, dims = dims_use, verbose = FALSE)
brain <- FindClusters(brain, resolution = 0.5, random.seed = 42, verbose = FALSE)
brain <- RunUMAP(brain, dims = dims_use, seed.use = 42, verbose = FALSE)

pdf(file.path(out_dir, "05_umap_clusters.pdf"), width = 7, height = 6)
print(DimPlot(brain, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()

# -----------------------------
# 7. Marker discovery
# -----------------------------
markers <- FindAllMarkers(
  brain,
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25
)

write.csv(markers, file.path(out_dir, "cluster_markers.csv"), row.names = FALSE)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

write.csv(top_markers, file.path(out_dir, "top10_markers_per_cluster.csv"), row.names = FALSE)

# -----------------------------
# 8. Broad annotation helper
# -----------------------------
marker_panels <- list(
  Excitatory_neuron = c("Slc17a7", "Rbfox3", "Snap25"),
  Inhibitory_neuron = c("Gad1", "Gad2", "Slc32a1"),
  Astrocyte = c("Gfap", "Aqp4", "Aldoc"),
  Oligodendrocyte = c("Mbp", "Mog", "Plp1"),
  OPC = c("Pdgfra", "Cspg4"),
  Microglia = c("Tmem119", "P2ry12", "Cx3cr1"),
  Endothelial = c("Kdr", "Pecam1", "Cldn5")
)

available <- unique(unlist(lapply(marker_panels, function(x) x[x %in% rownames(brain)])))

if (length(available) > 0) {
  pdf(file.path(out_dir, "06_marker_dotplot.pdf"), width = 11, height = 6)
  print(DotPlot(brain, features = available) + RotatedAxis())
  dev.off()
}

brain$cell_type_reviewed <- paste0("Cluster_", Idents(brain))

write.csv(
  brain@meta.data,
  file.path(out_dir, "cell_metadata.csv"),
  row.names = TRUE
)

saveRDS(brain, file.path(out_dir, "brain_seurat_object.rds"))

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

message("Done. Review marker tables before assigning final biological labels.")
