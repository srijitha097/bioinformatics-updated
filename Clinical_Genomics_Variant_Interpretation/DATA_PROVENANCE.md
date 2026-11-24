# Data provenance

The research workspace retains the original downloads and a phase-two cache under `comparative_genomics/data/`. Source file SHA-256 hashes are recorded in the [phase-two manifest](comparative_genomics/source_hashes.txt). Baseline file hashes are available in the [baseline manifest](comparative_genomics/baseline_hashes.txt).

| Phase-two input | Source and selection |
|---|---|
| Published phenopackets | [Phenopacket Store archive](https://codeload.github.com/monarch-initiative/phenopacket-store/tar.gz/refs/heads/main), retrieved 24 November 2026. Each selected archive member and its PMID are retained in [cases.tsv](comparative_genomics/rare/cases.tsv). Main is mutable. |
| Phenotype annotations | Baseline HPO gene-to-phenotype snapshot, hash in the rare-disease summary. For the publication-restricted audit, [HPO disease annotations](https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa), release 2026-09-02. Disease IDs link the annotation tables. |
| Population frequency | [gnomAD GraphQL](https://gnomad.broadinstitute.org/api), gnomad_r4. Thirty GRCh38 variant queries have archived individual responses. Fourteen lack available counts. Unknown frequency is retained separately from observed zero. |
| Germline expert assertions | [NCBI EFetch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi), db=clinvar, rettype=vcv, is_variationid=true. Two XML responses retain expert narratives and versioned SCVs. Selection intentionally sampled reference classes and is excluded from ranking validation. OMIM full text was unavailable. Disease IDs come through public annotation cross-references. |
| Somatic contexts and literature | CIViC 01-Nov-2026 accepted-only exports from phase one. Each [context](comparative_genomics/oncology/contexts.tsv) lists its AID, evidence ID, and PMID. PubMed EFetch, db=pubmed and retmode=xml, returned 23 publication records. One lacked an abstract. FDA sources are linked in the report and regulatory addendum. |
| Exome benchmark | [Zenodo v2](https://doi.org/10.5281/zenodo.3636192). Same GRCh37 NA12878 ARUP VCF, assay capture, matched GIAB v3.3.2 truth/confidence and OMIM gene spans as phase one. No new specimen pairing. |
| Allele reference | [1000 Genomes GRCh37 FASTA](https://s3.amazonaws.com/1000genomes/technical/reference/human_g1k_v37.fasta.gz). Indexed expanded FASTA with samtools 1.25. Independently tested selected autosomal allele records with bcftools norm 1.25. |
| Real CNV comparison | Baseline downsampled HCC1143 tumor/normal ICE slices. [CCLE GRCh37 segment file](https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/cellline_ccle_broad/data_cna_hg19.seg), HCC1143_BREAST, chromosome 20. Targets must fall wholly inside a segment. Cell line identity matches, while passage identity is unknown. |
| Simulated CNVs | Generated locally with development seed 504 and test seed 8841. No patient or array observations enter simulator development or evaluation. |

The normalized comparison requires the two phase-one `verify_*.bed` autosomal BEDs. The first-phase inputs, thresholds and reference links are described in each original project README.

## Redistribution and attribution

The NA12878 coordinates and derived calls use Erinija Pranckeviciene’s *NA12878 WES Benchmark dataset*, version 2. The panel files intersect or subset the deposited intervals. The record specifies CC BY 4.0. [Dataset and author](https://zenodo.org/records/3636192) · [License](https://creativecommons.org/licenses/by/4.0/)

The selected case table extracts identifiers and observed terms from the Monarch Initiative Phenopacket Store. Its BSD notice accompanies the extraction. CIViC releases its knowledgebase content under CC0. The literature tables identify the publications and summarize evidence in original prose. Raw publication records, submitted ClinVar narratives, HPO releases, gnomAD API responses, BAMs and reference files stay in the local input cache. [Phenopacket notice](comparative_genomics/rare/PHENOPACKET_LICENSE.txt) · [CIViC licensing](https://docs.civicdb.org/en/latest/about/faq.html)

The legacy CCLE segment mirror lacks a verified file-specific redistribution license. The public summary reports target and segment counts. Running the array comparison locally reconstructs the row-level join from the cited input. [CCLE data](https://sites.broadinstitute.org/ccle/datasets) · [DepMap licensing discussion](https://forum.depmap.org/t/re-distribution-of-the-ccle-data-licensing/3063)

To reconstruct the comparative cache, save the linked Phenopacket archive as `packets.tgz`, disease annotations as `phenotype.hpoa`, the legacy segment file as `ccle.seg`, and the reference as `b37.fa.gz`. Decompress the reference to `b37.fa` and run `samtools faidx comparative_genomics/data/b37.fa`. Use ClinVar EFetch with `db=clinvar&rettype=vcv&is_variationid=true&from_esearch=true` and the following comma-separated variation IDs. Save the first response as `clinvar.xml` and the second as `extra.xml`.

```text
17000,17023,7107,7109,17670,17676,17681,37389,42537,164113,164114,14093,42644,42744,164070,177780,161303,42325
576,15239
```

For `abstracts.xml`, use PubMed EFetch with `retmode=xml` and the PMIDs in the context table. The replay fetches uncached gnomAD responses by allele and builds `af.json`. Current database responses may differ from the archived XML or JSON. The supplied hash manifest identifies the exact inputs used here. [PMIDs](comparative_genomics/oncology/contexts.tsv) · [Archived input hashes](comparative_genomics/source_hashes.txt)
