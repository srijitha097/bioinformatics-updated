suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# Load the Seurat object created by the main analysis
brain <- readRDS("results/brain_seurat_object.rds")

# First-pass cell-type annotation based on cluster marker genes
celltype_map <- c(
  "0"  = "Excitatory neuron",
  "1"  = "SATB2+ excitatory neuron",
  "2"  = "Neural progenitor",
  "3"  = "NPY+ inhibitory neuron",
  "4"  = "Neuron",
  "5"  = "SST interneuron",
  "6"  = "Developing neuron",
  "7"  = "Developing inhibitory neuron",
  "8"  = "Striatal neuron",
  "9"  = "Neuron",
  "10" = "Intermediate neural progenitor",
  "11" = "CHODL+ inhibitory neuron",
  "12" = "HTR3A+ interneuron",
  "13" = "Cajal-Retzius-like neuron",
  "14" = "OPC",
  "15" = "Microglia"
)

brain <- RenameIdents(
  brain,
  celltype_map
)

# Save cell-type annotation in metadata
brain$cell_type <- Idents(brain)

# Plot annotated UMAP
p <- DimPlot(
  brain,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Cell types in embryonic mouse brain") +
  NoLegend()

ggsave(
  "results/07_umap_celltypes.pdf",
  p,
  width = 10,
  height = 7
)

# Save annotated metadata
write.csv(
  brain@meta.data,
  "results/cell_metadata_annotated.csv"
)

# Save annotated Seurat object
saveRDS(
  brain,
  "results/brain_seurat_annotated.rds"
)

cat("Cell-type annotation complete.\n")
cat("Annotated cell types:\n")
print(table(brain$cell_type))
