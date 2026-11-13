# HG002 Germline Variant Calling and Quality Assessment

## Overview

This analysis examines germline variation in a chromosome 20 region of the HG002 human genome using publicly available whole-genome sequencing data. I used an existing aligned Illumina sequencing dataset to assess read coverage, call variants with GATK HaplotypeCaller, normalize and filter VCF records, and compare the resulting alleles with Genome in a Bottle (GIAB) reference calls.

The analysis covers `chr20:10000000-11000000` (GRCh38; 1-based inclusive coordinates), approximately 1 Mb of the original whole-genome sequencing dataset.

## Data

- **Sample:** HG002 (NA24385), Illumina PCR-free whole-genome sequencing data, approximately 35× coverage.
- **Input alignment:** Preprocessed chromosome 20 BAM from the public DeepVariant case-study dataset.
- **Reference:** GRCh38/hg38 chromosome 20 FASTA.
- **Comparison data:** GIAB HG002 v4.2.1 GRCh38 benchmark VCF and high-confidence regions BED.

The aligned reads are used as input; alignment from original raw FASTQ files is not part of the main analysis.

## Methods

1. **Alignment QC:** Summarized the input BAM with SAMtools `flagstat` and calculated per-position depth across the target region.
2. **Germline variant calling:** Called SNPs and indels using GATK HaplotypeCaller.
3. **VCF processing:** Normalized alleles and split multiallelic records using BCFtools. Applied an exploratory `QUAL < 30` filter flag and summarized the resulting variants.
4. **Reference comparison:** Compared normalized REF/ALT alleles with GIAB calls within the benchmark high-confidence regions, including PASS or unfiltered calls.
5. **Reporting:** Generated a coverage summary, allele-overlap statistics, and an HTML QC report.

## Results

### Sequencing coverage

| Metric | Result |
|---|---:|
| Positions assessed | 1,000,001 |
| Mean depth | 34.994× |
| Positions with ≥10× coverage | 99.937% |
| Positions with ≥20× coverage | 99.187% |
| Positions with zero coverage | 0.008% |

### GIAB allele comparison

The high-confidence BED covered 977,840 bases within the selected region.

| Metric | Result |
|---|---:|
| Called alleles | 1,607 |
| GIAB reference alleles | 1,606 |
| Matching alleles (TP) | 1,605 |
| Unmatched called alleles (FP) | 2 |
| Missed reference alleles (FN) | 1 |
| Allele-level precision | 99.876% |
| Allele-level recall | 99.938% |
| Allele-level F1 score | 99.907% |

The called alleles closely matched the GIAB alleles in this region. These values describe **exact normalized allele overlap**, not genotype-level correctness or haplotype-aware benchmarking.

## Reproducing the analysis

Run the following commands from the repository root after installing the required tools (`samtools`, `bcftools`, GATK 4, Python 3, and Java):

```bash
bash scripts/check_tools.sh
bash scripts/get_data.sh
bash scripts/run.sh
python3 -m unittest discover -s tests -v
```

The data script first attempts to retrieve an indexed regional BAM. If the remote host cannot provide the required byte ranges, a larger chromosome 20 BAM download may be necessary:

```bash
FULL_BAM=1 bash scripts/get_data.sh
```

After the run, open the report and inspect the summary files:

```bash
open results/report.html
cat results/coverage_summary.tsv
cat results/benchmark.json
```

## Output files

| File | Contents |
|---|---|
| `results/flagstat.txt` | Alignment statistics |
| `results/coverage_summary.tsv` | Depth and coverage metrics |
| `results/variants.raw.vcf.gz` | Germline variants from HaplotypeCaller |
| `results/variants.norm.vcf.gz` | Normalized variant records |
| `results/variants.filtered.vcf.gz` | Variants with exploratory low-QUAL flags |
| `results/variant_stats.txt` | BCFtools variant summary |
| `results/benchmark.json` | GIAB allele comparison |
| `results/report.html` | HTML QC and comparison report |
| `results/versions.txt` | Recorded software versions |

The repository can include selected reports in `docs/results/`; large input files and intermediate sequencing data are excluded from version control.

## Interpretation and limitations

This was a regional analysis of a WGS-derived dataset, not variant calling across the full genome. The input BAM was already aligned and duplicate-processed. The reference comparison uses exact normalized REF/ALT matches within GIAB high-confidence regions and does not evaluate genotypes or complex haplotypes. Different representations of indels can affect matching, and the `QUAL < 30` filter is exploratory rather than a validated clinical filtering strategy.

## Data sources

- [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)
- [HG002 GIAB GRCh38 v4.2.1 benchmark files](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/)
- [UCSC hg38 reference sequences](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/)
