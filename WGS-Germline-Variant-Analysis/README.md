# Project 1 — Real HG002 WGS Region: Germline Calling, QC & Benchmarking

**Original engineering question:** Can a workflow take a real public human Illumina aligned dataset through sample-level QC, regional GATK HaplotypeCaller, illustrative QUAL-based filtering, variant normalization, and honest allele-level QC comparison with GIAB?

## Real data provenance

- Aligned, duplicate-processed Illumina 35x WGS data: Google DeepVariant case-study public **HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam**, from PrecisionFDA Truth Challenge (Google hosts chr20 only).
- UCSC hg38/GRCh38 chromosome 20 reference `chr20.fa.gz`; corresponding chr20 sequence must match BAM reference `M5`/reference hash. Check readme notes before benchmarking.
- GIAB HG002 GRCh38 v4.2.1 truth VCF and BED, NCBI: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/
- Analysis target `chr20:10000000-11000000`. Coordinates use 1-based inclusive regions in samtools/GATK.

## Mac prerequisites

A conda or micromamba environment with: `samtools`, `bcftools`, `gatk4`, `python >=3.10`, optionally `nextflow`, `bwa-mem2`, `fastqc`, `multiqc`. Tool availability depends on Mac architecture; conda Linux container may be easier on Apple Silicon. Examples:

```bash
# If you have conda with bioconda configured
conda create -n iu-wgs -c conda-forge -c bioconda \
  python=3.11 samtools bcftools gatk4 -y
conda activate iu-wgs
samtools --version | head -1
bcftools --version | head -1
gatk --version
```

## Run step by step

```bash
cd 01_wgs_germline_benchmark
bash scripts/check_tools.sh

# Downloads / remote region extraction. Several GB may be needed for the fallback.
# Samtools first tries indexed HTTP range extraction; if that fails, a full chr20 BAM
# download will be required. See printed fallback command.
bash scripts/get_data.sh

# Independent QC + variant calling + normalized allele overlap
bash scripts/run.sh

# Result outputs:
# results/flagstat.txt, coverage_summary.tsv, variants.raw.vcf.gz,
# variants.norm.vcf.gz, variants.filtered.vcf.gz, truth.norm.vcf.gz, benchmark.json, report.html

# Python logic tests (no downloaded data required)
python3 -m unittest discover -s tests -v
```

If remote indexed BAM access fails, download the complete chr20 BAM and matching `.bai` into `data/`, then re-run `get_data.sh` with `FULL_BAM=1` (will subset locally). See script for exact URLs. The chr20 BAM is **not** a 1 Mb BAM; it can be several gigabytes.

### Optional additional alignment exercise

After `get_data.sh`, optionally run `bash scripts/realign_optional.sh` (requires `bwa-mem2`, FastQC recommended). This reconstructs paired reads from the regional aligned BAM, runs FastQC, re-aligns them, and marks duplicates. Then run `BAM=data/HG002.realigned.dedup.bam bash scripts/run.sh`. The extracted FASTQs are **not original raw WGS reads**, and region boundary pairs may be lost. The default path keeps the documented preprocessed BAM for reproducibility.

### Nextflow option

Once you have data downloaded:

```bash
nextflow run main.nf --bam data/HG002.region.bam \
  --reference data/chr20.fa --truth data/HG002.truth.vcf.gz \
  --bed data/HG002.high_confidence.bed --region chr20:10000000-11000000
```

The Nextflow driver delegates to the same standard shell analysis and requires tools available on the execution host; it does not implicitly install them.

### Benchmarking caveat

`compare_vcfs.py` evaluates **exact REF/ALT allele matches after left-normalization** only in GIAB high-confidence regions; it does not resolve complex haplotypes and does not compare genotype correctness. A simplistic QUAL<30 exclusion is used as an exploratory flag only, not a clinically validated filtering procedure. Its precision/recall are **illustrative allele-level concordance**, not publication-grade benchmarking. For a more rigorous extension, use hap.py or RTG `vcfeval` with the complete matching reference.

### What to inspect

- `coverage_summary.tsv`: what proportion of positions have >=10 and >=20 read depth?
- `results/benchmark.json`: numerator/denominator counts and limitations.
- `results/report.html`: review warning thresholds.
- `results/versions.txt`: capture the program versions for reproducibility.

**Resume after successful run:** "Implemented an HG002 chr20 germline-variant QC workflow using samtools, GATK and BCFtools; automated allele-level GIAB comparison and HTML reporting with documented validation limitations." (Do not call the regional run whole-genome execution.)
