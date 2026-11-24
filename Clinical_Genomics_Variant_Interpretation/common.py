import csv
import gzip
import hashlib
import json
from pathlib import Path


def op(path):
    return (
        gzip.open(path, "rt", encoding="utf-8")
        if str(path).endswith(".gz")
        else open(path, encoding="utf-8")
    )


def rows(path):
    with op(path) as f:
        lines = (line for line in f if not line.startswith("##"))
        reader = csv.DictReader(lines, delimiter="\t")
        for row in reader:
            yield {
                k.lstrip("#").strip(): (v or "").strip()
                for k, v in row.items()
                if k is not None
            }


def write(path, records, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(records)


def report(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()
