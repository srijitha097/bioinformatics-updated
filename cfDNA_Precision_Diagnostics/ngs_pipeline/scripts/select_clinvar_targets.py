from pathlib import Path

import pandas as pd

DATA = Path("ngs_pipeline/data/clinvar")

df = pd.read_csv(
    DATA / "target_gene_variants.tsv",
    sep="\t",
    dtype=str
)

df["review_rank"] = 0

df.loc[
    df.review_status.str.contains(
        "multiple_submitters",
        na=False
    )
    & df.review_status.str.contains(
        "no_conflicts",
        na=False
    ),
    "review_rank"
] = 2

df.loc[
    df.review_status.str.contains(
        "reviewed_by_expert_panel",
        na=False
    ),
    "review_rank"
] = 3

df.loc[
    df.review_status.str.contains(
        "practice_guideline",
        na=False
    ),
    "review_rank"
] = 4

high_confidence = df[
    (df.clinical_significance == "Pathogenic")
    & (df.review_rank >= 2)
].copy()

high_confidence = high_confidence.sort_values(
    ["matched_gene", "review_rank", "position"],
    ascending=[True, False, True]
)

high_confidence.to_csv(
    DATA / "high_confidence_variants.tsv",
    sep="\t",
    index=False
)

print("\nHIGH-CONFIDENCE PATHOGENIC SNVs")

for gene in ["HBB", "CFTR", "HEXA", "GJB2"]:

    group = high_confidence[
        high_confidence.matched_gene
        .str.split("|")
        .apply(lambda genes: gene in genes)
    ]

    print(f"\n===== {gene} =====")
    print("High-confidence variants:", len(group))

    print(
        group[[
            "chromosome",
            "position",
            "clinvar_id",
            "ref",
            "alt",
            "rsid",
            "review_rank",
            "condition"
        ]].head(8).to_string(index=False)
    )

print("\nSPECIFIC TARGET VARIANTS")

targets = {
    "rs334",
    "rs75527207",
    "rs80338950"
}

selected = df[
    df.rsid.str.split("|").apply(
        lambda values: bool(
            set(values) & targets
        )
    )
]

print(selected.to_string(index=False))

print(
    "\nTotal high-confidence variants:",
    len(high_confidence)
)
