# Admixed Genotype Simulation
A Snakemake pipeline for simulating genotypes for admixed population. The aim is to have fully automatic pipeline for simulating genotypes for admixed populations.

## Note for updated software
**Please see https://kangchenghou.github.io/admix-kit/simulate-admix-genotype.html for a more intuitive and easy-to-use interface.**

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
chr_i = 2
kg_legend = pd.read_csv(f'/u/project/pasaniuc/pasaniucdata/admixture/1000G_haplotype/1000GP_Phase3/1000GP_Phase3_chr{chr_i}.legend.gz', delim_whitespace=True)
ukb_bim = pd.read_csv(f"/u/project/sgss/UKBB/data/cal/{chr_i}.bim", delim_whitespace=True, header=None)
ukb_index = kg_legend['position'].isin(ukb_bim[3].values)
duplicated = kg_legend['position'].duplicated()
biallelic = kg_legend['TYPE'] == "Biallelic_SNP"
np.savetxt(f"data/ukb_array.snp_list", kg_legend.loc[ukb_index & (~duplicated) & biallelic, 'id'].values, fmt='%s')
```

2. Extract a subsample of SNPs from 1kg
```python
import numpy as np
import pandas as pd
chr_i = 2
n_snp = 1_000
kg_legend = pd.read_csv(f'/u/project/pasaniuc/pasaniucdata/admixture/1000G_haplotype/1000GP_Phase3/1000GP_Phase3_chr{chr_i}.legend.gz', delim_whitespace=True)
maf_threshold = 0.01

maf_mode = "AND"
pops = ["AFR", "EUR"]
# filter biallelic SNPs and SNPs with the given MAF threshold in ALL / ANY population
maf_filter_index = ((maf_threshold < kg_legend[pops]) &
                    (kg_legend[pops] < 1 - maf_threshold))
if maf_mode == 'AND':
    maf_filter_index = maf_filter_index.all(axis=1)
elif maf_mode == 'OR':
    maf_filter_index = maf_filter_index.any(axis=1)

biallelic = kg_legend['TYPE'] == "Biallelic_SNP"
filter_index = np.where(maf_filter_index & biallelic)[0]
filter_index = filter_index[np.linspace(0, len(filter_index) - 1, num=n_snp + 2)[1:-1].astype(int)]

np.savetxt(f"data/kg_1k.snp_list", kg_legend.loc[filter_index, 'id'].values, fmt='%s')
```


1. Extract five regions of dense SNPs 
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
  
