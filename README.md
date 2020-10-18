# Admixed Genotype Simulation
A Snakemake pipeline for simulating genotypes for admixed population. The aim is to have fully automatic pipeline for simulating genotypes for admixed populations.

## Requirements
1. Reference phased genotype panels: download from https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html
2. admix-simu software: download from https://github.com/williamslab/admix-simu
3. HAPGEN2: download from https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html

## Steps
1. Snakemake rule `extract_raw`: extract ancestry reference population and perform basic quality control. We extract the ancestral population according to 1KG sample files. And filter for bi-allelic SNPs that have minor allele frequency larger than the specificied one in  <ins>all</ins> of the ancestral populations.
2. Snakemake rule `extend_hapgen`: extend the ancestral population using HAPGEN2.
3. Snakemake rule `simulate_admixture`: Simulate admixture population from the extended ancestral population.

## Run
First set up the configuration, specified in `config.yaml` properly. Then submit the script
```bash
qsub snakemake.sh
```

## To-do
- Improve documentation.
- Performance improvement at `extend_hapgen` and `simulate_admixture` steps.

  
## Improvement
The pipeline is in active delveopment, please post issue / pull request for suggestions.