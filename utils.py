import fire
import numpy as np
import pandas as pd
import subprocess
from os.path import join

import os
from shutil import copyfile
from itertools import compress
import gzip

def extract_raw_data(raw_dir, out_dir, pops, chr_i, maf_threshold, maf_mode):
    """
    extract_raw_data: extract data from 1000G Phase3 phased data,
        downloaded from https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html
        
        raw_dir: Raw data directory
        out_dir: Extracted data directory
        pops: which population to extract
        chr_i: chromosome to extract TODO: None to indicate all of the chromosomes
    """
    sample = pd.read_csv(join(raw_dir, '1000GP_Phase3.sample'), delim_whitespace=True)
    map_file = join(raw_dir, f'genetic_map_chr{chr_i}_combined_b37.txt')
    legend = pd.read_csv(join(raw_dir, f'1000GP_Phase3_chr{chr_i}.legend.gz'), delim_whitespace=True)
    with gzip.open(join(raw_dir, f'1000GP_Phase3_chr{chr_i}.hap.gz'), 'rt', encoding='utf-8') as f:
        haps = f.readlines()

    # filter biallelic SNPs and SNPs with the given MAF threshold in ALL / ANY population
    maf_filter_index = ((maf_threshold < legend[pops]) & 
                        (legend[pops] < 1 - maf_threshold))
    if maf_mode == 'AND':
        maf_filter_index = maf_filter_index.all(axis=1)
    elif maf_mode == 'OR':
        maf_filter_index = maf_filter_index.any(axis=1)
    filter_index = maf_filter_index & (legend['TYPE'] == 'Biallelic_SNP')
    print(f'#extract_raw_data: filtering, {sum(filter_index)} / {len(filter_index)} left')
    legend = legend[filter_index].reset_index(drop=True)
    haps = list(compress(haps, filter_index))
    haps = [hap.strip().replace(' ', '') for hap in haps]

    copyfile(map_file, join(out_dir, f'chr{chr_i}.map'))
    legend.to_csv(join(out_dir, f'chr{chr_i}.legend'), index=False, sep=' ')
    
    for pop in pops:

        pop_index = np.repeat((sample['GROUP'] == pop).values, 2)
        pop_haps = [' '.join(compress(hap, pop_index)) for hap in haps]
        with open(join(out_dir, f'chr{chr_i}.{pop}.hap'), 'w') as f:
            f.writelines('\n'.join(pop_haps))

    # produce SNP file here
    snp_info = {'SNP': legend['id'].values, 
            'CHR': chr_i,
            'M': legend['position'].values / 1e8,
            'POS': legend['position'].values,
            'A1': legend['a0'].values,
            'A2': legend['a1'].values}
    pd.DataFrame(snp_info).to_csv(join(out_dir, f'chr{chr_i}.snp'), header=False, index=False, sep=' ', float_format='%.6f')

def extend_hapgen(map_file, legend_file, hap_file, num_haplos, hapgen2_path, out_prefix):
    """
    extend_hapgen: extend the population
        map_file: path to the PLINK files for the seed population
        legend_file: path to the legend file
        hap_file: path to the haplotype file
        num_haplos: the desired number of haplotypes to generate
        hapgen_path: path of hapgen2
        out_prefix: prefix for the output
    """
    
    # cope with disease loci
    positionOfDiseaseSNP = pd.read_csv(legend_file, delim_whitespace=True)['position'][0]
    riskAllele = 1
    hetRisk = 1
    homRisk = 1

    # format the file properly for hapgen2 output
    subprocess.run([hapgen2_path, 
                '-m', map_file,
                '-l', legend_file,
                '-h', hap_file,
                '-o', out_prefix,
                '-dl', str(positionOfDiseaseSNP), str(riskAllele), str(hetRisk), str(homRisk),
                '-n', str(num_haplos), str(0),
                '-no_gens_output'])
    
    
    with open(f'{out_prefix}.controls.haps') as f:
        out_hap = [line.strip().replace(' ', '') for line in f.readlines()]

    with open(f'{out_prefix}.phgeno', 'w') as f:
        f.writelines('\n'.join(out_hap))
    

def simulate_admixture(input_prefix_list, admix_prop, num_gens, num_haplos, snp_file, out_dir, admix_simu_dir):
    """
    simulate_admixture:
        input_prefix_list: a #populations-length list
        admix_prop: a #populations-length list, relative proportion of the admixture, 
            will be rescaled to 1.
        num_gens: number of generations to simulate
        num_haplos: number of haplotypes to simulate
        out_prefix: output prefix
        admix_simu_dir: Directory to the admix_simu software
    """
    
    assert (len(input_prefix_list) == len(admix_prop)), "`input_prefix_list`, `admix_prop` should have the same length"
    num_pops = len(input_prefix_list)
    assert (num_pops == 2), "Currently we only support two-way admixture"

    # read the SNP file, confirm that the two snp file are the same
    phgeno_list = [f'-POP{i + 1} {input_prefix}.phgeno' for i, input_prefix in enumerate(input_prefix_list)]


    admix_dat = ['\t'.join([str(num_haplos), 'ADMIX', *[f'POP{i}' for i in np.arange(1, num_pops + 1)]]),
                 '\t'.join([str(num_gens), '0', *[str(prop) for prop in admix_prop]])]

    dat_file = join(out_dir, 'admix.dat')
    with open(dat_file, 'w') as f:
        f.writelines('\n'.join(admix_dat))

    cmd = ' '.join([join(admix_simu_dir, 'simu-mix-2n.pl'),
                    dat_file,
                    snp_file,
                    join(out_dir, 'admix'),
                    *phgeno_list])
    print(cmd)
    os.system(cmd)

    cmd = ' '.join([join(admix_simu_dir, 'bp2anc.pl'),
                    join(out_dir, 'admix.bp'),
                    '>',
                    join(out_dir, 'admix.hanc')])
    print(cmd)
    os.system(cmd)
