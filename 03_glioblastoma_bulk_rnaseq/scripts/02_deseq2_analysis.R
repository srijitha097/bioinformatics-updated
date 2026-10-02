suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(tibble)
})

set.seed(42)

dir.create("results", recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(
  "data/processed/raw_count_matrix.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.delim(
  "data/processed/sample_metadata.tsv",
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$gsm
metadata$condition <- factor(
  metadata$condition,
  levels = c("normal", "glioblastoma")
)

counts <- counts[, rownames(metadata), drop = FALSE]

# Remove very low-count genes.
keep <- rowSums(counts >= 10) >= 10
counts_filtered <- counts[keep, , drop = FALSE]

message("Genes before filtering: ", nrow(counts))
message("Genes after filtering:  ", nrow(counts_filtered))

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filtered),
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "glioblastoma", "normal"),
  alpha = 0.05
)

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

write.csv(
  res_df,
  "results/deseq2_all_genes.csv",
  row.names = FALSE
)

sig <- res_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )

write.csv(
  sig,
  "results/deseq2_significant_genes.csv",
  row.names = FALSE
)

message("Significant DE genes: ", nrow(sig))

# ---------------------------
# PCA
# ---------------------------
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(
  pca_data,
  aes(PC1, PC2, shape = condition, label = name)
) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(
  "results/pca_gbm_vs_normal.pdf",
  p_pca,
  width = 7,
  height = 5
)

# ---------------------------
# MA plot
# ---------------------------
pdf("results/ma_plot.pdf", width = 7, height = 5)
plotMA(res, alpha = 0.05, ylim = c(-6, 6))
dev.off()

# ---------------------------
# Volcano plot
# ---------------------------
volcano_df <- res_df %>%
  mutate(
    neglog10_padj = -log10(pmax(padj, 1e-300)),
    status = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >= 1 ~ "Up",
      !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

p_volcano <- ggplot(
  volcano_df,
  aes(log2FoldChange, neglog10_padj)
) +
  geom_point(
    aes(shape = status),
    alpha = 0.55,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    x = "log2 fold change: glioblastoma vs normal",
    y = "-log10 adjusted p-value",
    title = "Glioblastoma vs normal brain"
  ) +
  theme_bw()

ggsave(
  "results/volcano_plot.pdf",
  p_volcano,
  width = 7,
  height = 5
)

# ---------------------------
# Heatmap: top 40 DE genes
# ---------------------------
top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 40) %>%
  pull(gene_id)

mat <- assay(vsd)[top_genes, , drop = FALSE]
mat_z <- t(scale(t(mat)))

annotation_col <- data.frame(
  condition = metadata$condition
)
rownames(annotation_col) <- rownames(metadata)

pdf("results/top40_deg_heatmap.pdf", width = 11, height = 10)
pheatmap(
  mat_z,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "Top 40 differentially expressed genes"
)
dev.off()

# Save normalized counts for inspection.
norm_counts <- counts(dds, normalized = TRUE)
write.csv(
  as.data.frame(norm_counts),
  "results/deseq2_normalized_counts.csv"
)

capture.output(
  sessionInfo(),
  file = "results/sessionInfo.txt"
)
