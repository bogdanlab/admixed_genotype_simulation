from os.path import join
from utils import extract_raw_data, extend_hapgen, simulate_admixture
import os

configfile: "config.yaml"


rule all:
    input:
        join("out/0_raw/chr{}.snp".format(config['CHR'])),
        expand("out/1_hapgen/chr{}.{{pop}}.phgeno".format(config['CHR']), pop=config["POPS"]),
        "out/2_admix_geno/EUR_0.5_AFR_0.5_10_10000/admix.phgeno"

rule extract_raw:
    input:
        join(config["INPUT_DATA_DIR"], 'genetic_map_chr{}_combined_b37.txt'.format(config['CHR'])),
        join(config["INPUT_DATA_DIR"], '1000GP_Phase3_chr{}.legend.gz'.format(config['CHR'])),
        join(config["INPUT_DATA_DIR"], '1000GP_Phase3_chr{}.hap.gz'.format(config['CHR']))
    output:
        "out/0_raw/chr{}.snp".format(config['CHR']),
        "out/0_raw/chr{}.map".format(config['CHR']),
        "out/0_raw/chr{}.legend".format(config['CHR']),
        expand("out/0_raw/chr{}.{{pop}}.hap".format(config['CHR']), pop=config["POPS"])
    run:
        if not os.path.exists("out/0_raw"):
            os.makedirs("out/0_raw")
        extract_raw_data(raw_dir=config["INPUT_DATA_DIR"], 
                         out_dir="out/0_raw", 
                         pops=config['POPS'],
                         chr_i=config['CHR'],
                         maf_threshold=config['MAF_THRESHOLD'])

rule extend_hapgen:
    input:
        map_file="out/0_raw/chr{}.map".format(config['CHR']),
        snp_file="out/0_raw/chr{}.snp".format(config['CHR']),
        legend_file="out/0_raw/chr{}.legend".format(config['CHR']),
        hap_file="out/0_raw/chr{}.{{pop}}.hap".format(config['CHR'])
    output:
        "out/1_hapgen/chr{}.{{pop}}.phgeno".format(config['CHR'])
    run:
        if not os.path.exists("out/1_hapgen"):
            os.makedirs("out/1_hapgen")
        
        extend_hapgen(map_file=input.map_file,
                      legend_file=input.legend_file,
                      hap_file=input.hap_file,
                      num_haplos=config['HAPGEN_NUM_HAPLOS'],
                      hapgen2_path=config['HAPGEN2_PATH'],
                      out_prefix=output[0][0 : output[0].rfind('.')])
        

rule simulate_admixture:   
    input:
        snp_file="out/0_raw/chr{}.snp".format(config['CHR']),
        pop1_phgeno="out/1_hapgen/chr{}.{{pop1}}.phgeno".format(config['CHR']),
        pop2_phgeno="out/1_hapgen/chr{}.{{pop2}}.phgeno".format(config['CHR']),
    output:
        "out/2_admix_geno/{pop1}_{prop1}_{pop2}_{prop2}_{num_gens}_{num_haplos}/admix.phgeno"
    run:
        out_dir = output[0][0 : output[0].rfind('/')]
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        pop1_prefix = input.pop1_phgeno.rsplit('.', 1)[0]
        pop2_prefix = input.pop2_phgeno.rsplit('.', 1)[0]
        
        if wildcards.num_haplos >= config['HAPGEN_NUM_HAPLOS'] / 2:
            print("""
                  Warning: number of haplotypes in the admixture population should be smaller than 
                  half of ancestral population haplotypes. Please adjust numbers accordingly.
                  """)
        
        simulate_admixture(input_prefix_list=[pop1_prefix, pop2_prefix], 
                           admix_prop=[wildcards.prop1, wildcards.prop2], 
                           num_gens=wildcards.num_gens,
                           num_haplos=wildcards.num_haplos, 
                           snp_file=input.snp_file,
                           out_dir=out_dir, 
                           admix_simu_dir=config['ADMIX_SIMU_DIR'])
        

