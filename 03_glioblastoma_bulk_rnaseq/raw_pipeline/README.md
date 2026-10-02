# Optional raw FASTQ workflow

The main repository begins with the public GEO raw count files because this makes the analysis practical for code review.

The source study also provides raw FASTQs in SRA.

A full upstream workflow would be:

```text
SRA FASTQ
   |
FastQC / MultiQC
   |
HISAT2 -> GRCh38
   |
SAMtools sort/index
   |
HTSeq-count
   |
raw gene-count matrix
   |
DESeq2
```

The original GEO study reports using HISAT2 v2.1.0 with GRCh38 and HTSeq-count v0.10.0.

A production pipeline should pin the genome/annotation release and software versions in its environment or container.
