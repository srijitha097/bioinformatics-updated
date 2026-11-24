import json
import tempfile
import unittest
from pathlib import Path

from common import rows, write
from exome_gene_panel_cnv_analysis.run import compare
from oncology_somatic_variant_interpretation.run import audit, extract
from rare_disease_variant_curation.run import prioritize


class LiveSchemaRegressions(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def test_hpo_does_not_merge_terms_across_diseases(self):
        p = self.root
        variant = dict(
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
        )
        write(p / "variants.tsv", [variant], list(variant))
        hpo = [
            dict(gene_symbol="CFTR", hpo_id=term, disease_id=disease)
            for term, disease in [("HP:0002110", "OMIM:1"), ("HP:0001738", "OMIM:2")]
        ]
        hpo = [dict(row, frequency="-") for row in hpo]
        hpo.append(dict(hpo[1], disease_id="OMIM:1", frequency="HP:0040285"))
        write(p / "hpo.tsv", hpo, list(hpo[0]))
        write(
            p / "cases.tsv",
            [{"case_id": "fixture", "hpo_ids": "HP:0002110;HP:0001738"}],
            ["case_id", "hpo_ids"],
        )
        prioritize(p / "variants.tsv", p / "hpo.tsv", p / "cases.tsv", p, 1)
        candidate = next(rows(p / "candidate_queue.tsv"))
        self.assertEqual(float(candidate["hpo_overlap"]), 0.5)
        self.assertEqual(candidate["hpo_best_disease_ids"], "OMIM:1;OMIM:2")

    def test_civic_requires_explicit_statusless_export_and_matching_pubmed(self):
        p = self.root
        assertion = dict(
            assertion_id="1",
            molecular_profile="Fixture",
            disease="Fixture",
            assertion_type="Diagnostic",
            amp_category="Tier II - Level C",
            evidence_item_ids="5",
        )
        write(p / "assertions.tsv", [assertion], list(assertion))
        with self.assertRaisesRegex(ValueError, "No accepted"):
            extract(p / "assertions.tsv", p)
        extract(p / "assertions.tsv", p, accepted_only=True)
        evidence = dict(
            evidence_id="5",
            evidence_status="accepted",
            source_type="PubMed",
            citation_id="123",
        )
        write(p / "evidence.tsv", [evidence], list(evidence))
        review = next(rows(p / "review_template.tsv"))
        review.update(
            reviewed_pmids="456",
            evidence_summary="Fixture",
            reviewer="Fixture",
            review_date="2026-11-24",
            research_tier="Tier II - Level C",
            reviewed_evidence_ids="5",
        )
        write(p / "review.tsv", [review], list(review))
        with self.assertRaisesRegex(ValueError, "PMIDs"):
            audit(
                p / "accepted_assertions.tsv", p / "review.tsv", p, p / "evidence.tsv"
            )
        review["reviewed_pmids"] = "123"
        write(p / "review.tsv", [review], list(review))
        audit(p / "accepted_assertions.tsv", p / "review.tsv", p, p / "evidence.tsv")

    def test_benchmark_contig_aliases_and_half_open_boundary(self):
        p = self.root
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        (p / "truth.vcf").write_text(header + "1\t3\t.\tA\tG\t60\tPASS\t.\tGT\t0/1\n")
        (p / "query.vcf").write_text(
            header + "chr1\t3\t.\tA\tG\t60\tPASS\t.\tGT\t1|0\n"
            "chr1\t4\t.\tA\tG\t60\tPASS\t.\tGT\t0/1\n"
        )
        (p / "target.bed").write_text("chr1\t0\t3\nchr1\t1\t3\n")
        (p / "confidence.bed").write_text("1\t2\t5\n")
        compare(
            p / "truth.vcf",
            p / "query.vcf",
            p / "target.bed",
            p / "confidence.bed",
            p / "result.json",
        )
        result = json.loads((p / "result.json").read_text())
        self.assertEqual(
            (result["target_bases"], result["confident_target_bases"]), (3, 1)
        )
        self.assertEqual(
            (
                result["true_positive"],
                result["false_positive"],
                result["false_negative"],
            ),
            (1, 0, 0),
        )
