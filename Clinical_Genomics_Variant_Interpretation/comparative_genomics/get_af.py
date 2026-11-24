import csv
import json
import subprocess

from ranking import P

af = {}
for row in csv.DictReader(open(P / "rare/cases.tsv"), delimiter="\t"):
    allele = row["allele"]
    path = P / "data" / f"af_{allele}.json"
    if not path.exists():
        query = (
            "{ variant(variantId: "
            + json.dumps(allele)
            + ", dataset: gnomad_r4) { exome { ac an } genome { ac an } } }"
        )
        subprocess.run(
            [
                "curl",
                "-fsSL",
                "--retry",
                "2",
                "--max-time",
                "25",
                "https://gnomad.broadinstitute.org/api",
                "-H",
                "Content-Type: application/json",
                "--data",
                json.dumps({"query": query}),
                "-o",
                str(path),
            ],
            check=True,
        )
    response = json.loads(path.read_text())
    for error in response.get("errors", []):
        if error.get("message") != "Variant not found":
            raise ValueError(f"gnomAD query failed for {allele}: {error}")
    d = response.get("data") or {}
    counts = [
        v["ac"] / v["an"] for v in (d.get("variant") or {}).values() if v and v["an"]
    ]
    af[allele] = max(counts, default=None)
(P / "data/af.json").write_text(json.dumps(af, indent=2))
