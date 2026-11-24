"""Encode literature judgments before joining reference tiers. See selection caveats in README."""

import csv
import xml.etree.ElementTree as ET

from ranking import P, save


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# Assessments concern evidence available in this limited review.
ONCOLOGY = {
    5: (
        "defer",
        "The cited paper emphasizes uncommon EGFR alterations. Its abstract does not establish L858R-specific erlotinib benefit.",
    ),
    8: (
        "II-D",
        "One metastatic sarcoma with LMNA::NTRK1 regressed on LOXO-101. This case alone cannot estimate response across partner genes and tumor types.",
    ),
    10: (
        "I-A",
        "The phase 1b combination study suggested activity in BRAF V600 melanoma. The US cobimetinib label explicitly includes V600E with vemurafenib in unresectable or metastatic melanoma.",
    ),
    11: (
        "I-A",
        "The randomized coBRIM trial supports the combination in V600-mutant advanced melanoma. The FDA indication explicitly covers V600K, although the pooled estimate is not allele-specific.",
    ),
    20: (
        "II-C",
        "Stage IV colorectal cancers with V600E had lower observed two-year survival than D594G and wild-type groups. A retrospective comparison is susceptible to treatment and disease-subtype differences.",
    ),
    21: (
        "II-C",
        "Fusion FISH supported identification of mesenchymal chondrosarcoma in adequately prepared specimens. Inadequate samples lacked signal, limiting negative inference.",
    ),
    24: (
        "II-C",
        "A recurrent chimeric transcript was detected in 15 of 15 fibrolamellar tumors. A discovery cohort does not establish diagnostic specificity across all liver tumors.",
    ),
    33: (
        "II-D",
        "Gilteritinib produced antileukemic activity in FLT3-mutant model systems. Model response alone gives limited clinical evidence for the D835 subgroup.",
    ),
    34: (
        "II-C",
        "A single-arm phase 1/2 trial observed alectinib responses in ALK-rearranged advanced NSCLC. Lack of a randomized comparator limits this paper alone.",
    ),
    38: (
        "I-B",
        "ADMIRAL found longer survival versus salvage chemotherapy in relapsed/refractory FLT3-mutated AML. The trial pooled FLT3 alterations, limiting inference about ITD-specific effect size.",
    ),
    43: (
        "defer",
        "The Ph-like ALL study establishes a class of kinase lesions. The abstract lacks sufficient SNX2::ABL1-specific diagnostic discrimination.",
    ),
    51: (
        "defer",
        "The source title reports a pediatric NTRK2-fused tumor response. No abstract was available to confirm the exact partner or response assessment.",
    ),
    52: (
        "defer",
        "Two phase I trials showed activity in fusion-selected cancers. The abstract does not resolve the clinical effect for EML4::NTRK3 specifically.",
    ),
    56: (
        "II-D",
        "One infantile fibrosarcoma with the named fusion responded radiographically to LOXO-101. A single observation cannot quantify durability across the disease.",
    ),
    76: (
        "II-C",
        "Tumor breakpoint and transcript studies support a recurrent EWSR1::WT1 association. Application as a diagnostic marker requires integration with morphology.",
    ),
    79: (
        "II-C",
        "CIC::DUX4 characterized a subset of EWSR1-negative primitive sarcomas. Selection of the study cohort limits specificity estimates.",
    ),
    80: (
        "II-C",
        "BCOR internal tandem duplications recurred in childhood clear cell sarcoma of kidney. The molecular association supports a diagnostic workup within this histology.",
    ),
    85: (
        "I-A",
        "The phase II study found little activity with KIT D816V. US imatinib labeling limits the aggressive systemic mastocytosis indication to absent or unknown D816V.",
    ),
    87: (
        "II-C",
        "A series of non-RELA ZFTA fusions supports rearrangement testing in supratentorial ependymomas. Generalization to every possible fusion partner is uncertain.",
    ),
    93: (
        "defer",
        "Methylation-based reclassification identified CNS neuroblastoma with FOXR2 activation. The abstract does not validate rearrangement alone as a diagnostic test.",
    ),
    97: (
        "II-C",
        "A round-cell tumor study compared FISH and RT-PCR. Transcript testing can identify the partner, whereas break-apart testing by itself cannot identify FLI1.",
    ),
    98: (
        "II-C",
        "Integrated testing of suspected Ewing sarcoma identified rearrangements and fusion partners. RNA quality and primers affect detection of individual subtypes.",
    ),
    99: (
        "II-D",
        "One subcutaneous tumor carrying the fusion showed atypical immunostaining. This observation supports a differential diagnosis with substantial uncertainty.",
    ),
}


def oncology():
    papers = {
        a.findtext("PMID"): a
        for a in ET.parse(P / "data/abstracts.xml").findall(".//MedlineCitation")
    }
    out = []
    for r in read(P / "oncology/contexts.tsv"):
        tier, note = ONCOLOGY[int(r["aid"])]
        title = "".join(papers[r["pmid"]].find("Article/ArticleTitle").itertext())
        out.append(dict(r, assessment=tier, reason=note, publication=title))
    save("oncology/assessments.tsv", out)
    ref = {
        a["assertion_id"]: a["amp_category"]
        for a in read(
            next(
                Path("oncology_somatic_variant_interpretation/data").glob(
                    "*AcceptedAssertion*"
                )
            )
        )
    }
    save("oncology/comparison.tsv", [dict(r, reference=ref[r["aid"]]) for r in out])
    (P / "oncology/report.md").write_text(
        "# Molecular oncology evidence reviews\n\nResearch only. Tier and level describe the evidence examined. Selection used accepted CIViC contexts and linked studies. Some reference tiers had already been visible during phase one. Judgments were recorded separately before the comparison join.\n\n"
        + "\n\n".join(
            f"## {r['profile']} in {r['disease']}\n\n{r['context']} · {r['therapy'] or r['significance']} · {r['assessment']}. {r['reason']} [Primary publication](https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/) · [CIViC context](https://civicdb.org/links/assertions/{r['aid']})"
            for r in out
        )
        + "\n\n[US cobimetinib label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=c387579e-cee0-4334-bd1e-73f93ac1bde6) · [US imatinib label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=211ef2da-2868-4a77-8055-1cb2cd78e24b)\n"
    )


from pathlib import Path

GERM = {
    "576": (
        "PVS1, PS3, PM3, PP4_Moderate",
        "The PAH panel cited a canonical splice donor, functional loss, recessive observations and phenotype specificity. Functional evidence may overlap with the splice mechanism; this reconstructs panel reasoning.",
    ),
    "15239": (
        "PP1_Strong, PS4_Moderate, PM3, PM2_Supporting, PP3, PS3_Supporting",
        "The hemoglobinopathy panel used segregation, case enrichment, trans observations and rarity. It gave the splice assay supporting strength.",
    ),
    "7107": (
        "CFTR2 historical disease-liability assessment",
        "The truncating allele was classified by CFTR2 for cystic fibrosis. The submission links PMID 23974870 but has no explicit ACMG codes. Individual-level functional data were unavailable here.",
    ),
    "17000": (
        "PS4, PP1_Strong, PM3, PP3",
        "Two strong criteria support pathogenicity in the historical hearing-loss submission. Reduced penetrance complicates frequency cutoffs. Earlier baseline used another panel record with different strengths.",
    ),
    "17023": (
        "PS4, PP1_Strong, PM3",
        "Case enrichment and strong segregation meet a pathogenic combination despite the high population frequency.",
    ),
    "17670": (
        "ENIGMA multifactorial",
        "Very low posterior probability and European population frequency support the historical class 1 assessment.",
    ),
    "17676": (
        "ENIGMA multifactorial",
        "Panel posterior probability 2.23e-6 supports class 1. Avoid counting its component observations again.",
    ),
    "17681": (
        "ENIGMA multifactorial",
        "Panel posterior probability below 1e-10 supports class 1 under the stated IARC thresholds.",
    ),
    "37389": (
        "ENIGMA multifactorial",
        "Low missense and splicing priors together with ancestry-specific population frequency supported the panel interpretation.",
    ),
    "42537": (
        "PS4_Moderate, PM2_Supporting, PP1_Strong, PM1",
        "Segregation with case enrichment and regional evidence supports likely pathogenicity for dominant HCM.",
    ),
    "164113": (
        "PS4_Moderate, PP1_Strong, PM1",
        "Segregation and recurrence support likely pathogenicity. The panel did not require computational support.",
    ),
    "164114": (
        "PS4_Moderate, PM2_Supporting, PM1, PP3",
        "Case enrichment and the affected region contribute moderate evidence, supplemented by rarity and prediction.",
    ),
    "14093": (
        "PS4_Moderate, PM6, PM2, PP3",
        "Multiple moderate criteria support likely pathogenicity under the historical MYH7 specification.",
    ),
    "42644": (
        "PVS1_Moderate, PS4_Supporting, PP1_Supporting",
        "The splice position alone does not justify full PVS1. The retained strengths leave uncertainty.",
    ),
    "42744": (
        "PS4_Moderate, PM2_Supporting, PVS1_Moderate",
        "A terminal truncation limits loss-of-function strength. Tissue observations did not justify PS3.",
    ),
    "164070": (
        "PS4_Supporting, PP3",
        "Case observations and prediction are insufficient for likely pathogenicity.",
    ),
    "177780": (
        "BP7",
        "Synonymous status and limited splice evidence do not supply a second independent benign criterion.",
    ),
    "161303": ("BS1", "The VCEP permits BS1 alone for likely benign in this setting."),
    "42325": (
        "BS1, BP2, BP5",
        "Population frequency and the observed clinical contexts support the FBN1 panel assessment for Marfan syndrome.",
    ),
}


def germline():
    reviews = []
    for v in [
        *ET.parse(P / "data/clinvar.xml").getroot(),
        *ET.parse(P / "data/extra.xml").getroot(),
    ]:
        vid = v.get("VariationID")
        if vid not in GERM:
            continue
        a = next(
            a
            for a in v.findall(".//ClinicalAssertion")
            if a.findtext("Classification/ReviewStatus") == "reviewed by expert panel"
            and a.findtext("Classification/GermlineClassification", "").lower()
            in {
                "pathogenic",
                "benign",
                "uncertain significance",
                "likely pathogenic",
                "likely benign",
            }
        )
        scv = a.find("ClinVarAccession")
        codes, note = GERM[vid]
        reviews.append(
            dict(
                vid=vid,
                hgvs=v.get("VariationName"),
                scv=scv.get("Accession") + "." + scv.get("Version", ""),
                condition=" | ".join(
                    [
                        x.text or ""
                        for x in a.findall("TraitSet/Trait/Name/ElementValue")
                    ]
                    + [
                        x.get("DB", "") + ":" + x.get("ID", "")
                        for x in a.findall("TraitSet/Trait/XRef")
                    ]
                ),
                criteria=codes,
                reconstruction=a.findtext("Classification/GermlineClassification"),
                interpretation=note,
                pmids=";".join(
                    sorted(
                        {
                            x.text
                            for x in a.findall(".//Citation/ID")
                            if x.get("Source") == "PubMed" and x.text
                        }
                    )
                ),
                source="https://www.ncbi.nlm.nih.gov/clinvar/variation/" + vid + "/",
            )
        )
    save("rare/evidence_reviews.tsv", reviews)


if __name__ == "__main__":
    germline()
    oncology()
