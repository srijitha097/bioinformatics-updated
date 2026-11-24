# Laboratory QC and Assay Monitoring

Centralized sequencing quality assessment, automated failure detection, and interactive assay performance monitoring for cell-free DNA sequencing workflows.

## Monitoring Features

- Sample-level sequencing quality assessment.
- Mapping efficiency and mapping quality monitoring.
- Fragment length distribution analysis.
- Sequencing read count consistency.
- Configurable quality thresholds.
- Automated technical failure detection.
- Molecular consensus quality monitoring.
- Fetal fraction warnings.
- Posterior confidence and no-call monitoring.
- Interactive visualization and downloadable reports.

## Alignment QC Summary

| Metric | Result |
|---|---:|
| Assays monitored | 2 |
| Sequencing samples | 215 |
| Sequencing reads | 23,000,332 |
| Samples passing alignment QC | 215 |
| Samples failing alignment QC | 0 |

## Automated Failure Detection

Controlled quality failures were introduced into temporary copies of the sequencing benchmark data.

| Failure Scenario | Validation |
|---|---|
| Low mapping efficiency | Passed |
| Low pairing rate | Passed |
| Low high-quality mapping rate | Passed |
| Low mean mapping quality | Passed |
| Short fragment length | Passed |
| Long fragment length | Passed |
| Read count mismatch | Passed |
| Upstream QC failure | Passed |
| Missing mapping information | Passed |

All nine injected failure scenarios were detected without false alerts among the 206 unchanged samples.

## Assay Monitoring

| Event Category | Events |
|---|---:|
| Molecular benchmark warnings | 7 |
| Low fetal fraction warnings | 43 |
| Single-gene no-calls | 662 |
| Aneuploidy no-calls | 12 |
| Technical failure events | 0 |
| Total monitoring events | 724 |

Technical failures, benchmark warnings, and informational no-call events are classified separately.

## Generate Monitoring Reports

    python laboratory_qc/scripts/aggregate_qc.py
    python laboratory_qc/scripts/assay_monitoring.py

## Validate Automated Failure Detection

    python laboratory_qc/scripts/test_qc_alerts.py

## Launch Dashboard

    python -m streamlit run laboratory_qc/dashboard.py

Open http://localhost:8501 in a web browser.

## Implementation

- [QC threshold configuration](config/thresholds.json)
- [Centralized QC aggregation](scripts/aggregate_qc.py)
- [Failure detection validation](scripts/test_qc_alerts.py)
- [Assay-specific monitoring](scripts/assay_monitoring.py)
- [Interactive dashboard](dashboard.py)
- [Monitoring results](../results/laboratory_qc/)

## Technologies

Python, Pandas, Streamlit, Plotly, JSON, and automated QC validation.

## Research Scope

Monitoring thresholds are experimental. Ground-truth-dependent molecular metrics are used for synthetic benchmarking rather than clinical laboratory decision-making.
