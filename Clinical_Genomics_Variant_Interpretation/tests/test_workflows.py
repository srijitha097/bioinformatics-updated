import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from common import rows, write
from exome_gene_panel_cnv_analysis.cnv_demo import demo
from exome_gene_panel_cnv_analysis.run import compare
from oncology_somatic_variant_interpretation.run import audit as oncology_audit
from oncology_somatic_variant_interpretation.run import extract
from rare_disease_variant_curation.run import audit, prioritize


class Tests(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.TemporaryDirectory()
        self.p = Path(self.t.name)

    def tearDown(self):
        self.t.cleanup()

    def test_rare(self):
        p = self.p
        v = dict(
            Assembly="GRCh38",
            Type="single nucleotide variant",
            GeneSymbol="CFTR",
            VariationID="1",
            ReviewStatus="reviewed by expert panel",
            ClinicalSignificance="Pathogenic",
            ReferenceAlleleVCF="A",
            AlternateAlleleVCF="G",
            PositionVCF="10",
            Chromosome="7",
            Name="synthetic variant",
            PhenotypeList="synthetic disease",
        )
        write(p / "v.tsv", [v, dict(v, VariationID="2", Assembly="GRCh37")], list(v))
        write(
            p / "h.tsv",
            [dict(gene_symbol="CFTR", hpo_id="HP:0002110")],
            ["gene_symbol", "hpo_id"],
        )
        write(
            p / "c.tsv",
            [dict(case_id="x", hpo_ids="HP:0002110")],
            ["case_id", "hpo_ids"],
        )
        prioritize(p / "v.tsv", p / "h.tsv", p / "c.tsv", p, 5)
        r = list(rows(p / "candidate_queue.tsv"))
        self.assertEqual(len(r), 1)
        self.assertEqual(float(r[0]["hpo_overlap"]), 1)
        q = list(rows(p / "curation_template.tsv"))
        q[0].update({k: "synthetic fixture only" for k in q[0] if not q[0][k]})
        q[0]["research_classification"] = "Pathogenic"
        write(p / "review.tsv", q, list(q[0]))
        audit(p / "candidate_queue.tsv", p / "review.tsv", p)
        self.assertEqual(
            json.loads((p / "curation_summary.json").read_text())["reviewed"], 1
        )

    def test_oncology(self):
        p = self.p
        r = dict(
            assertion_id="1",
            assertion_status="accepted",
            amp_asco_cap_category="Tier I",
            disease="synthetic tumor",
            molecular_profile="synthetic variant",
            assertion_type="Predictive",
        )
        write(
            p / "a.tsv",
            [r, dict(r, assertion_id="2", assertion_status="submitted")],
            list(r),
        )
        extract(p / "a.tsv", p)
        self.assertEqual(len(list(rows(p / "accepted_assertions.tsv"))), 1)
        q = list(rows(p / "review_template.tsv"))
        q[0].update({k: "synthetic fixture only" for k in q[0] if not q[0][k]})
        q[0]["research_tier"] = "Tier I"
        write(p / "review.tsv", q, list(q[0]))
        oncology_audit(p / "accepted_assertions.tsv", p / "review.tsv", p)
        self.assertEqual(
            json.loads((p / "review_summary.json").read_text())["reviewed"], 1
        )

    def test_snv(self):
        p = self.p
        header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        line = lambda n, gt: f"chr1\t{n}\t.\tA\tG\t60\tPASS\t.\tGT\t{gt}\n"
        (p / "t.vcf").write_text(header + line(5, "0/1") + line(7, "1/1"))
        (p / "q.vcf").write_text(
            header + line(5, "1|0") + line(7, "0/1") + line(30, "0/1")
        )
        (p / "b.bed").write_text("chr1\t0\t10\n")
        compare(p / "t.vcf", p / "q.vcf", p / "b.bed", p / "b.bed", p / "o.json")
        d = json.loads((p / "o.json").read_text())
        self.assertEqual(
            (d["true_positive"], d["false_positive"], d["false_negative"]), (1, 1, 1)
        )
        self.assertEqual(d["f1"], 0.5)

    def test_coverage(self):
        p = self.p
        (p / "b.bed").write_text("chr1\t0\t4\nchr1\t2\t5\n")
        (p / "d.tsv").write_text("chr1\t1\t20\nchr1\t2\t30\n")
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "exome_gene_panel_cnv_analysis/coverage.py"),
                "--bed",
                str(p / "b.bed"),
                "--depth",
                str(p / "d.tsv"),
                "--out",
                str(p / "o.json"),
            ],
            check=True,
        )
        d = json.loads((p / "o.json").read_text())
        self.assertEqual(d["target_bases"], 5)
        self.assertEqual(d["mean_depth"], 10)
        self.assertEqual(d["fraction_20x"], 0.4)

    def test_cnv(self):
        demo(self.p)
        d = json.loads((self.p / "cnv_summary.json").read_text())
        self.assertTrue(d["synthetic_only"])
        self.assertEqual(len(list(rows(self.p / "synthetic_targets.tsv"))), 1200)
        self.assertGreater(d["recall"], 0.5)


if __name__ == "__main__":
    unittest.main()
