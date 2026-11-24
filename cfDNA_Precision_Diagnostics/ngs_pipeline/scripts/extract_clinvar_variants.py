import csv
import re
from collections import Counter
from pathlib import Path

import pysam
import pandas as pd

DATA = Path("ngs_pipeline/data/clinvar")

SOURCE = DATA / "clinvar.vcf.gz"
OUTPUT = DATA / "target_gene_variants.tsv"
VCF_OUTPUT = DATA / "target_gene_variants.vcf.gz"

TARGET_GENES = {
    "HBB",
    "CFTR",
    "HEXA",
    "GJB2"
}

VALID_SIGNIFICANCE = {
    "pathogenic",
    "likely_pathogenic",
    "pathogenic/likely_pathogenic"
}


def info_text(record, key):

    value = record.info.get(key)

    if value is None:
        return "."

    if isinstance(value, (tuple, list)):
        return "|".join(str(v) for v in value)

    return str(value)


rows = []
gene_counts = Counter()
classification_counts = Counter()
excluded = Counter()

with pysam.VariantFile(str(SOURCE)) as vcf:

    with pysam.VariantFile(
        str(VCF_OUTPUT),
        "wz",
        header=vcf.header
    ) as out_vcf:

        for record in vcf:

            gene_info = info_text(
                record, "GENEINFO"
            )

            genes = {
                item.split(":")[0].strip()
                for item in re.split(r"[|,]", gene_info)
                if ":" in item
            }

            matched_genes = genes & TARGET_GENES

            if not matched_genes:
                continue

            for gene in matched_genes:
                gene_counts[gene] += 1

            significance = info_text(
                record, "CLNSIG"
            )

            normalized = (
                significance.lower()
                .replace(" ", "_")
            )

            classification_counts[significance] += 1

            if normalized not in VALID_SIGNIFICANCE:
                excluded["clinical_significance"] += 1
                continue

            if (
                record.alts is None
                or len(record.alts) != 1
            ):
                excluded["not_biallelic"] += 1
                continue

            ref = record.ref.upper()
            alt = record.alts[0].upper()

            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in "ACGT"
                or alt not in "ACGT"
                or ref == alt
            ):
                excluded["not_snv"] += 1
                continue

            rsid = info_text(record, "RS")

            if rsid != ".":
                rsid = "|".join(
                    value
                    if value.startswith("rs")
                    else "rs" + value
                    for value in re.split(r"[|,]", rsid)
                )

            rows.append({
                "chromosome": record.chrom,
                "position": record.pos,
                "clinvar_id": record.id,
                "ref": ref,
                "alt": alt,
                "gene_info": gene_info,
                "clinical_significance": significance,
                "review_status": info_text(
                    record, "CLNREVSTAT"
                ),
                "condition": info_text(
                    record, "CLNDN"
                ),
                "rsid": rsid,
                "matched_gene": "|".join(
                    sorted(matched_genes)
                )
            })

            out_vcf.write(record)

if not rows:
    raise RuntimeError(
        "No variants retained. "
        "Inspect the annotation summaries."
    )

df = pd.DataFrame(rows)

df.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)

print("\nTOTAL GENE-ASSOCIATED VARIANTS")

for gene in sorted(TARGET_GENES):
    print(gene, gene_counts[gene])

print("\nCLINICAL SIGNIFICANCE VALUES")

for classification, count in (
    classification_counts.most_common(20)
):
    print(classification, count)

print("\nEXCLUDED VARIANTS")

for reason, count in excluded.most_common():
    print(reason, count)

print("\nTOTAL RETAINED VARIANTS:", len(df))

for gene in sorted(TARGET_GENES):

    group = df[
        df.matched_gene.str.split("|")
        .apply(lambda values: gene in values)
    ]

    print(f"\n===== {gene} =====")
    print("Retained pathogenic SNVs:", len(group))

    print("\nClinical significance:")

    print(
        group.clinical_significance
        .value_counts()
        .to_string()
    )

    print("\nExample variants:")

    print(
        group[[
            "chromosome",
            "position",
            "clinvar_id",
            "ref",
            "alt",
            "rsid",
            "clinical_significance",
            "review_status"
        ]].head(5).to_string(index=False)
    )

print("\nOutput:", OUTPUT)
print("Filtered VCF:", VCF_OUTPUT)
