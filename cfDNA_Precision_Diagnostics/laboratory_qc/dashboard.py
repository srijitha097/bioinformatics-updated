import json
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]

RESULTS = ROOT / "results"
QC = RESULTS / "laboratory_qc"
NGS = RESULTS / "ngs"

st.set_page_config(
    page_title="cfDNA Laboratory QC",
    page_icon="🧬",
    layout="wide"
)

st.title("cfDNA Laboratory QC & Assay Monitoring")

st.caption(
    "Research-only monitoring of synthetic "
    "prenatal sequencing experiments."
)


@st.cache_data
def load_csv(path):
    return pd.read_csv(path)


def load_json(path):
    with open(path) as handle:
        return json.load(handle)


qc = load_csv(
    QC / "centralized_alignment_qc.csv"
)

alerts = load_csv(
    QC / "qc_alerts.csv"
)

st.sidebar.header("Monitoring Filters")

assays = st.sidebar.multiselect(
    "Select assay",
    sorted(qc.assay.unique()),
    default=sorted(qc.assay.unique())
)

filtered = qc[
    qc.assay.isin(assays)
]

if filtered.empty:
    st.warning("Select at least one assay.")
    st.stop()

passed = (
    filtered.monitoring_status == "PASS"
).sum()

failed = (
    filtered.monitoring_status == "FAIL"
).sum()

cols = st.columns(5)

cols[0].metric(
    "Samples",
    f"{len(filtered):,}"
)

cols[1].metric(
    "Sequencing reads",
    f"{filtered.total_primary_reads.sum():,}"
)

cols[2].metric(
    "QC passed",
    int(passed)
)

cols[3].metric(
    "QC failed",
    int(failed)
)

cols[4].metric(
    "Pass rate",
    f"{passed / len(filtered):.2%}"
)

tabs = st.tabs([
    "Sequencing QC",
    "Molecular Quantification",
    "Fetal Fraction",
    "Aneuploidy Screening",
    "Single-Gene Screening",
    "QC Alerts"
])

with tabs[0]:

    st.subheader("Sequencing Quality Assessment")

    c1, c2 = st.columns(2)

    with c1:

        fig = px.box(
            filtered,
            x="assay",
            y="mapping_pct",
            color="assay",
            points="all",
            title="Mapping Efficiency"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    with c2:

        fig = px.box(
            filtered,
            x="assay",
            y="mean_fragment_length",
            color="assay",
            title="cfDNA Fragment Length"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    st.dataframe(
        filtered,
        use_container_width=True
    )

    st.download_button(
        "Download sample QC report",
        filtered.to_csv(index=False),
        file_name="sample_qc_report.csv",
        mime="text/csv"
    )


with tabs[1]:

    st.subheader("UMI Molecular Quantification")

    path = (
        NGS / "single_gene/benchmark/"
        "molecular_quantification_metrics.json"
    )

    m = load_json(path)

    cols = st.columns(4)

    cols[0].metric(
        "Original molecules",
        f"{m['original_molecules']:,}"
    )

    cols[1].metric(
        "Accepted families",
        f"{m['accepted_ref_alt_families']:,}"
    )

    cols[2].metric(
        "UMI families merged",
        f"{m['merged_umi_families']:,}"
    )

    cols[3].metric(
        "VAF error reduction",
        f"{m['vaf_error_reduction_pct']:.2f}%"
    )

    d = load_csv(
        NGS / "single_gene/benchmark/"
        "molecular_performance_by_depth.csv"
    )

    long = d.melt(
        id_vars=["gene", "molecules"],
        value_vars=[
            "raw_vaf_mae",
            "consensus_vaf_mae"
        ],
        var_name="method",
        value_name="MAE"
    )

    long["MAE_pp"] = long.MAE * 100

    fig = px.bar(
        long,
        x="gene",
        y="MAE_pp",
        color="method",
        barmode="group",
        facet_col="molecules",
        title="Raw versus Consensus VAF Error",
        labels={"MAE_pp": "MAE (percentage points)"}
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

    st.dataframe(d, use_container_width=True)


with tabs[2]:

    st.subheader("Fetal Fraction Estimation")

    path = (
        NGS / "fetal_fraction/"
        "fetal_fraction_summary.csv"
    )

    if path.exists():

        ff = load_csv(path)

        st.dataframe(
            ff,
            use_container_width=True
        )

        image = (
            NGS / "fetal_fraction/"
            "fetal_fraction_benchmark.png"
        )

        if image.exists():
            st.image(str(image))

    else:

        st.warning(
            "Fetal fraction summary unavailable."
        )


with tabs[3]:

    st.subheader("Chromosomal Aneuploidy Screening")

    m = load_json(
        NGS / "aneuploidy/benchmark/"
        "screening_metrics.json"
    )

    raw = m["raw"]
    corrected = m["gc_corrected"]

    cols = st.columns(4)

    cols[0].metric(
        "Called accuracy",
        f"{corrected['accuracy_called']:.2%}",
        f"{(corrected['accuracy_called'] - raw['accuracy_called'])*100:.2f} pp"
    )

    cols[1].metric(
        "Trisomy sensitivity",
        f"{corrected['exact_trisomy_sensitivity_called']:.2%}"
    )

    cols[2].metric(
        "Called specificity",
        f"{corrected['specificity_called']:.2%}"
    )

    cols[3].metric(
        "No-call rate",
        f"{corrected['no_call_rate']:.2%}"
    )

    d = load_csv(
        NGS / "aneuploidy/benchmark/"
        "euploid_baseline.csv"
    )

    long = d.melt(
        id_vars="chromosome",
        value_vars=[
            "raw_cv_pct",
            "corrected_cv_pct"
        ],
        var_name="method",
        value_name="CV"
    )

    fig = px.bar(
        long,
        x="chromosome",
        y="CV",
        color="method",
        barmode="group",
        title="GC Correction and Baseline Variability"
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

    calls = load_csv(
        NGS / "aneuploidy/benchmark/"
        "screening_calls.csv"
    )

    st.dataframe(
        calls,
        use_container_width=True
    )


with tabs[4]:

    st.subheader("Single-Gene Genotype Screening")

    m = load_json(
        NGS / "single_gene/benchmark/"
        "genotype_screening_metrics.json"
    )

    mode = st.radio(
        "Parental information",
        ["carrier_known", "paternal_unknown"],
        horizontal=True
    )

    metrics = m[mode]

    cols = st.columns(4)

    cols[0].metric(
        "Called accuracy",
        f"{metrics['accuracy_called']:.2%}"
    )

    cols[1].metric(
        "No-call rate",
        f"{metrics['no_call_rate']:.2%}"
    )

    cols[2].metric(
        "Affected detection",
        f"{metrics['affected_detection_rate']:.2%}"
    )

    cols[3].metric(
        "False positives",
        metrics["false_positives"]
    )

    d = load_csv(
        NGS / "single_gene/benchmark/"
        "genotype_performance_by_ff_depth.csv"
    )

    d = d[d["mode"] == mode].copy()

    d["Fetal fraction (%)"] = d.true_ff * 100

    d["Detection rate (%)"] = (
        d.affected_detection_rate * 100
    )

    d["Depth"] = d.input_depth.astype(str)

    fig = px.line(
        d.sort_values("true_ff"),
        x="Fetal fraction (%)",
        y="Detection rate (%)",
        color="Depth",
        markers=True,
        title="Affected Genotype Detection"
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

    st.dataframe(
        d,
        use_container_width=True
    )


with tabs[5]:

    st.subheader("Automated Laboratory QC Alerts")

    active = alerts[
        alerts.assay.isin(assays)
    ]

    if active.empty:

        st.success(
            "No QC threshold violations detected."
        )

    else:

        st.error(
            f"{len(active)} QC violations detected."
        )

        st.dataframe(
            active,
            use_container_width=True
        )

    st.download_button(
        "Download QC alert log",
        active.to_csv(index=False),
        file_name="qc_alerts.csv",
        mime="text/csv"
    )

st.divider()

st.caption(
    "Synthetic sequencing benchmarks using human "
    "reference sequences and ClinVar annotations. "
    "Not intended for clinical diagnostic use."
)
