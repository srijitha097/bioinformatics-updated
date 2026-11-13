# Publish Project 1 after analyzing real data

1. Run the real pipeline and check the report in `results/RUN_SUMMARY.md`.
2. Fill in `docs/RESULTS_TEMPLATE.md` with verified numbers, covariates, effect sizes, and limitations.
3. Cite NCBI GEO GSE254569 and the originating data-generation work, even though the analysis and code are original.
4. Keep raw and donor-linked measurements out of GitHub. `.gitignore` already excludes `data/raw`, `data/processed`, and `results`.
5. To publish reviewed, non-sensitive figures, copy only validated figures into a new `figures_public/` directory and clearly caption their origin.

```bash
cd ~/Documents/01_schizophrenia_snrnaseq
git init
git add README.md LICENSE environment.yml scripts tests docs run_demo.sh run_real.sh .gitignore
git status  # inspect every staged file first
git commit -m "Add reproducible human brain snRNA-seq workflow"
# Create an empty GitHub repo using GitHub.com, then follow its URL/remote instructions.
```

Never commit `donor_mapping_template.csv` if it contains private or non-public clinical metadata. Do not add metrics from synthetic tests to your resume.
