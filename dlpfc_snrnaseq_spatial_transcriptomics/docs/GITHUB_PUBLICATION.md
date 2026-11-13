# Publish Project 2 after analyzing real data

1. Run the real LIBD workflow and inspect `results/RUN_SUMMARY.md` and figures.
2. Complete `docs/RESULTS_TEMPLATE.md` and report sections/donors separately.
3. Cite official LIBD `spatialLIBD` resources and the underlying data-generation study (their references are linked in README).
4. Never publish cached single-nucleus or Visium raw donor matrices. The repo `.gitignore` excludes processed and raw data.
5. Place only reviewed, appropriately anonymized and non-sensitive figures in `figures_public/`, then explicitly label them real-data analysis.

```bash
cd ~/Documents/02_dlpfc_spatial_integration
git init
git add README.md LICENSE environment.yml scripts tests docs run_demo.sh run_real.sh .gitignore
git status
git commit -m "Add reproducible snRNA-seq and spatial transcriptomics workflow"
# Create an empty GitHub repo, then add its remote URL and push.
```
