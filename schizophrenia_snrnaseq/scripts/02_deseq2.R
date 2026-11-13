#!/usr/bin/env Rscript
## Pseudobulk DESeq2 with donor-celltype coverage and covariate checks.
## Intent: exploratory secondary analysis of public GSE254569 counts.
suppressPackageStartupMessages(library(DESeq2))
scriptarg <- grep('^--file=', commandArgs(trailingOnly=FALSE), value=TRUE)
if (!length(scriptarg)) stop('Run this via Rscript')
script <- sub('^--file=', '', scriptarg[1])
root <- normalizePath(file.path(dirname(script), '..'), mustWork=TRUE)
manifest <- read.csv(file.path(root, 'results', 'pseudobulk_manifest.csv'), stringsAsFactors=FALSE)
out <- file.path(root, 'results', 'deseq2')
dir.create(out, showWarnings=FALSE, recursive=TRUE)
min_nuclei <- 20L              # diagnostic QC threshold, not a universal statistical cutoff
min_donors_per_group <- 15L    # require adequate independent case/control replication
for (j in seq_len(nrow(manifest))) {
  ct <- manifest$folder[j]
  cat('Analyzing:', ct, '\n')
  path <- file.path(root, 'results', 'pseudobulk', ct)
  ctmat <- read.csv(file.path(path, 'counts.csv'), row.names=1, check.names=FALSE)
  md <- read.csv(file.path(path, 'metadata.csv'), stringsAsFactors=FALSE)
  if (anyDuplicated(md$sample_id)) stop('Duplicate pseudobulk sample IDs in ', ct)
  rownames(md) <- md$sample_id
  stopifnot(identical(colnames(ctmat), rownames(md)))
  if (!'n_nuclei' %in% names(md)) stop('Rerun upgraded prepare_pseudobulk: n_nuclei is missing')
  valid <- is.finite(md$n_nuclei) & md$n_nuclei >= min_nuclei & colSums(ctmat) > 0
  cat('  donor-celltype samples passing >=', min_nuclei, 'nuclei:', sum(valid), '/', nrow(md), '\n')
  md <- md[valid, , drop=FALSE]
  ctmat <- ctmat[, valid, drop=FALSE]
  md$diagnosis <- factor(md$diagnosis, levels=c('control','case'))
  nc <- sum(md$diagnosis=='case'); nn <- sum(md$diagnosis=='control')
  if (nc < min_donors_per_group || nn < min_donors_per_group) {
    warning(sprintf('Skipping %s: only %d case and %d control donors pass coverage QC', ct, nc, nn))
    writeLines(sprintf('SKIPPED %s: case donors=%d; control donors=%d; required minimum=%d; minimum nuclei/sample=%d',
      ct, nc, nn, min_donors_per_group, min_nuclei), file.path(out,paste0(ct,'_summary.txt')))
    next
  }
  keep <- rowSums(ctmat >= 10) >= 3
  ctmat <- as.matrix(ctmat[keep, , drop=FALSE])
  if (nrow(ctmat) < 10) {
    warning('Skipping ', ct, ': too few expressed genes')
    writeLines(paste('SKIPPED', ct, ': too few expressed genes'),file.path(out,paste0(ct,'_summary.txt')))
    next
  }
  storage.mode(ctmat) <- 'integer'
  candidates <- character()
  # Center and scale all continuous variables. Do not treat integer-valued age or PMI as categories.
  for (v in c('age','pmi','rin')) {
    if (!v %in% names(md)) next
    z <- suppressWarnings(as.numeric(md[[v]]))
    if (anyNA(z) || !all(is.finite(z)) || length(unique(z))<2) next
    vz <- paste0(v,'_z')
    md[[vz]] <- as.numeric(scale(z))
    candidates <- c(candidates,vz)
  }
  if ('sex' %in% names(md) && !anyNA(md$sex) && length(unique(md$sex))>1) {
    md$sex <- factor(md$sex); candidates <- c(candidates,'sex')
  }
  if ('lib_batch' %in% names(md) && !anyNA(md$lib_batch) && length(unique(md$lib_batch))>1) {
    md$lib_batch <- factor(md$lib_batch); candidates <- c(candidates,'lib_batch')
  }
  dvars <- c('diagnosis')
  for (v in candidates) {
    trial <- reformulate(c(dvars, v))
    mm <- model.matrix(trial, md)
    # Avoid rank deficient covariates and underdetermined model; batch is only included if estimable.
    if (qr(mm)$rank==ncol(mm) && nrow(mm)>=ncol(mm)+12L) dvars <- c(dvars,v)
    else cat('  dropping non-estimable covariate:',v,'\n')
  }
  formula <- reformulate(dvars)
  cat('  design:', deparse(formula), '\n')
  dds <- DESeqDataSetFromMatrix(countData=ctmat, colData=md, design=formula)
  dds <- DESeq(dds, quiet=TRUE)
  res <- results(dds, contrast=c('diagnosis','case','control'), alpha=0.05)
  tab <- data.frame(gene=rownames(res), as.data.frame(res), row.names=NULL)
  tab <- tab[order(tab$padj, na.last=TRUE), ]
  write.csv(tab, file.path(out,paste0(ct,'_deseq2.csv')), row.names=FALSE)
  png(file.path(out,paste0(ct,'_volcano.png')),width=1000,height=700,res=120)
  ok <- is.finite(tab$log2FoldChange) & is.finite(tab$padj)
  cols <- ifelse(tab$padj<0.05 & abs(tab$log2FoldChange)>=0.5, 'firebrick','gray70')
  plot(tab$log2FoldChange[ok], -log10(pmax(tab$padj[ok],1e-300)), pch=16, cex=.55,
    col=cols[ok], xlab='log2 fold change (case vs control)', ylab='-log10(FDR)',main=paste('Donor-level DE:',ct))
  abline(h=-log10(0.05),lty=2);dev.off()
  writeLines(paste('Cell type:',ct,'\nDonors:',ncol(ctmat),
    '\nCase donors after nucleus QC:',nc,'\nControl donors after nucleus QC:',nn,
    '\nMin nuclei/sample:',min_nuclei, '\nGenes tested:',nrow(ctmat),
    '\nDE genes (FDR < 0.05):',sum(tab$padj<0.05,na.rm=TRUE),
    '\nDesign:',deparse(formula)), file.path(out,paste0(ct,'_summary.txt')))
}
cat('Differential expression complete in results/deseq2/\n')
