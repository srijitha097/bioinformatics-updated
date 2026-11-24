# Progress

All steps below ran on 24 November 2026.

## Workspace and baseline

- Inspected the top-level and project READMEs, scripts, input fixtures and tests. The workspace contained three starter projects and no Git metadata. All five baseline tests passed.
- Checked hardware before downloads. The machine has 16 GiB of RAM and the workspace volume reported 109 GiB available. Retrieved the 459 MB compressed ClinVar summary after checking the size. Used existing local Conda binaries for plotting and sequencing tools without changing those environments.

## Project 1

- Retrieved HPO gene-to-phenotype annotations and streamed the GRCh38 ClinVar variant summary. Ran the original ranking as a baseline.
- Changed ranking to keep phenotype evidence within a single gene–disease annotation and excluded 3,083 zero-frequency HPO rows. Reran on the original fixed scenarios and added tests for cross-disease leakage and excluded annotations.
- Reconstructed two published ClinGen Hearing Loss VCEP assessments for GJB2 in DFNB1A. Verified both GRCh38 alleles in ClinVar, retrieved exome/genome gnomAD v4 allele counts, and cross-checked the OMIM disease identifier through the GeneReviews OMIM table. Direct OMIM access was unavailable.
- Ran the curation audit and checked HPO ties, the gene-level figure, and population count denominators.

## Project 2

- Reproduced the starter’s live CIViC extraction failure on the statusless combined export. Used CIViC GraphQL release metadata to verify the dated accepted-only assertion and evidence paths.
- Downloaded and parsed both 1 November 2026 files. Added an explicit mode for the statusless accepted-only export. Preserved explicit status filtering for files with status columns.
- Reviewed six selected assertions under the AMP/ASCO/CAP framework. Executed the audit with a join from accepted assertion to accepted evidence item and PubMed identifier. Added regression coverage for missing status provenance and mismatched PMIDs.
- Generated the source-tier table and plot. Checked that counts sum to the extraction totals.

## Project 3

- Downloaded the ARUP NA12878 exome callset, SeqCap EZ targets, GIAB v3.3.2 truth and confident regions, and OMIM gene spans. All five MD5 hashes matched Zenodo.
- Reproduced an empty truth intersection from UCSC BED versus GRCh37 VCF chromosome aliases. Verified the build, harmonized aliases, and computed callable benchmark space. Reran the autosomal exome and ten-gene subset benchmarks.
- Recomputed filtered genotype keys with bcftools 1.25 and confirmed identical TP, FP and FN counts for both scopes. Added a regression test covering the alias and BED boundary behavior.
- Downloaded the small matched HCC1143 tumor/normal exome slices and ICE intervals from the GATK public bucket. Checked BAM integrity, converted capture intervals from one-based inclusive to BED, and ran samtools depth. Measured coverage and generated exploratory relative copy ratios for targets meeting the depth floor.
- Replaced the idealized CNV demo’s known baseline with an estimate from separately simulated copy-neutral controls. Reran with held-out simulated cases and variable event widths. Inspected missed short events and false gains.
- Generated and reviewed the scientific figures, structured results, research reports and reproduction instructions. All eight tests passed.

## Second-phase research run · 24 November 2026

- Left first-phase results in place and created a separate comparative_genomics workspace and baseline checksum manifest. Rechecked disk space before retrieving the reference genome.
- Retrieved the Phenopacket Store archive, selected 30 solved cases, fetched actual gnomAD counts and ran ranking and perturbation experiments. Repeated phenotype scoring with HPO annotations backed by explicitly different publications. Recorded the marked loss of coverage and recovery.
- First ClinVar XML request returned empty records. Retried using explicit variation IDs, inspected individual expert SCVs and reconstructed 19 reviews. Some submissions lacked named conditions and used external condition IDs, which are preserved.
- Selected 23 new CIViC contexts and downloaded primary publication metadata. Completed restricted evidence assessments, joined reference tiers afterward, and recorded the selective post-comparison regulatory review separately. One source had no abstract.
- Corrected a syntax error in the phase-two benchmark driver. Verified the GRCh37 masks, normalization and independent exact-genotype counts. The FASTA gzip utility reported trailing bytes, so indexing was run separately and REF checks completed on the relevant capture records. The initial VCF header warnings persisted for unused legacy fields.
- Ran the expanded panel and depth/GQ sensitivities and wrote an allele-level discrepancy table. The expanded panel had zero qualifying truth indels, precluding recall measurement.
- Evaluated simulated CNV methods with separate development and test seeds. The first CCLE array download returned HTTP 404 and a study archive URL returned HTTP 403. Located the GRCh37 segments in the cBioPortal data repository and performed the matched-cell-line intersection.
- Generated three phase-two plots, expanded validation tests and checked baseline hashes. The projects remain research work in progress.
- Ran physical BAM subsampling at four fractions, using three seeds per reduced fraction. Full-input target coverage reproduced the saved slice measurements. The platform bytecode cache was unwritable, so syntax compilation was repeated successfully with a temporary cache.

## Final scientific audit

Replayed saved-input workflows and figures, verified all 15 tests, recalculated the restricted ranking metrics and reconciled oncology context discrepancies. Independently checked CNV true-target identities and locked-model reproducibility across seeds. Refreshed regulatory source links in the expanded reports. Final findings and claim-level evidence are indexed in [FINAL_AUDIT.md](FINAL_AUDIT.md).

## Repository quality review

Replayed the three studies in a temporary checkout with archived inputs. The first replay exposed samtools-dependent coverage differences. Pinning samtools and bcftools to 1.25 reproduced all baseline hashes. Rechecked the eight figures at the plotted-value level, added vector exports, and indexed the release files. [Review and verification](FINAL_GITHUB_REVIEW.md)
