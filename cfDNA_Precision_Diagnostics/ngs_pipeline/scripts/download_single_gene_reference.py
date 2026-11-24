import subprocess
import time
from pathlib import Path
from urllib.parse import urlencode

import pandas as pd

OUT = Path(
    "ngs_pipeline/data/reference/single_gene"
)

OUT.mkdir(parents=True, exist_ok=True)

targets = pd.read_csv(
    "ngs_pipeline/data/clinvar/selected_clinvar_targets.tsv",
    sep="\t"
)

for row in targets.itertuples():

    params = {
        "db": "nuccore",
        "id": row.refseq_accession,
        "seq_start": row.region_start_1based,
        "seq_stop": row.region_end_1based,
        "rettype": "fasta",
        "retmode": "text"
    }

    url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/efetch.fcgi?"
        + urlencode(params)
    )

    path = OUT / f"{row.gene}.fa"

    print(
        "Downloading", row.gene,
        row.refseq_accession,
        flush=True
    )

    subprocess.run([
        "curl",
        "-fL",
        "--retry", "5",
        "--connect-timeout", "30",
        url,
        "-o", str(path)
    ], check=True)

    subprocess.run([
        "samtools",
        "faidx",
        str(path)
    ], check=True)

    time.sleep(1)

print("\nAll reference regions downloaded")
