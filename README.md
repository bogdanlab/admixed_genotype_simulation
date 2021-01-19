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

## Helper scripts

1. Extract the list of SNPs from 1kg that's overlap with UKB
```python
import numpy as np
import pandas as pd
chr_i = 22
kg_legend = pd.read_csv(f'/u/project/pasaniuc/pasaniucdata/admixture/1000G_haplotype/1000GP_Phase3/1000GP_Phase3_chr{chr_i}.legend.gz', delim_whitespace=True)
ukb_bim = pd.read_csv(f"/u/project/sgss/UKBB/data/cal/{chr_i}.bim", delim_whitespace=True, header=None)
ukb_index = kg_legend['position'].isin(ukb_bim[3].values)
duplicated = kg_legend['position'].duplicated()
biallelic = kg_legend['TYPE'] == "Biallelic_SNP"
np.savetxt(f"data/ukb_array.{chr_i}.txt", kg_legend.loc[ukb_index & (~duplicated) & biallelic, 'id'].values, fmt='%s')
```

2. Extract five regions of SNPs 
```python
import numpy as np
import pandas as pd
chr_i = 22
kg_legend = pd.read_csv(f'/u/project/pasaniuc/pasaniucdata/admixture/1000G_haplotype/1000GP_Phase3/1000GP_Phase3_chr{chr_i}.legend.gz', delim_whitespace=True)
ukb_bim = pd.read_csv(f"/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/genotype/raw/chr{chr_i}.bim", delim_whitespace=True, header=None)
ukb_index = kg_legend['position'].isin(ukb_bim[3].values)
partition = pd.read_csv(f"/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/partition/fourier_ls-chr{chr_i}.bed", delim_whitespace=True)
block_start, block_stop = partition.loc[10, "start"], partition.loc[15, "stop"]
duplicated = kg_legend['position'].duplicated()
biallelic = kg_legend['TYPE'] == "Biallelic_SNP"
snp_list = kg_legend.loc[ukb_index & (~duplicated) & 
              biallelic & 
              (block_start <= kg_legend['position']) &
              (kg_legend['position'] < block_stop), 'id'].values
np.savetxt(f"data/ukb_imputed_blocks.txt", snp_list, fmt='%s')
```
  
