suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
})

dir.create("results", recursive = TRUE, showWarnings = FALSE)

de <- read.csv(
  "results/deseq2_significant_genes.csv",
  stringsAsFactors = FALSE
)

if (nrow(de) == 0) {
  stop("No significant DE genes found.")
}

# Remove Ensembl version suffix if present
de$ensembl <- sub("\\..*$", "", de$gene_id)

# Separate upregulated and downregulated genes
up_genes <- de %>%
  filter(log2FoldChange >= 1)

down_genes <- de %>%
  filter(log2FoldChange <= -1)

cat("Upregulated genes:", nrow(up_genes), "\n")
cat("Downregulated genes:", nrow(down_genes), "\n")

run_go <- function(gene_df, label) {

  mapped <- bitr(
    unique(gene_df$ensembl),
    fromType = "ENSEMBL",
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Hs.eg.db
  )

  write.csv(
    mapped,
    paste0("results/", label, "_gene_mapping.csv"),
    row.names = FALSE
  )

  ego <- enrichGO(
    gene = unique(mapped$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )

  go_df <- as.data.frame(ego)

  write.csv(
    go_df,
    paste0("results/GO_", label, ".csv"),
    row.names = FALSE
  )

  if (nrow(go_df) > 0) {
    pdf(
      paste0("results/GO_", label, "_dotplot.pdf"),
      width = 9,
      height = 7
    )

    print(
      dotplot(
        ego,
        showCategory = 20,
        title = paste("GO Biological Process -", label)
      )
    )

    dev.off()
  }

  cat(label, "GO terms:", nrow(go_df), "\n")
}

run_go(up_genes, "upregulated")
run_go(down_genes, "downregulated")
