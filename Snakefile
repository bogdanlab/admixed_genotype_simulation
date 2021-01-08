from os.path import join
from utils import *
import os
import yaml

configfile: "config.yaml"

with open("sim_setting.yaml", 'r') as f:
    sim_setting = yaml.safe_load(f)
    sim_prefix_list = [setting2prefix(*parse_sim_setting(sim)) for sim in sim_setting]

rule all:
    input:
        join("out/0_raw/chr{}.snp".format(config['CHR'])),
        expand("out/1_hapgen/chr{}.{{pop}}.phgeno".format(config['CHR']), pop=config["POPS"]),
        expand("out/2_admix_geno/{prefix}/admix.phgeno", prefix=sim_prefix_list)

rule extract_raw:
    input:
        join(config["INPUT_DATA_DIR"], 'genetic_map_chr{}_combined_b37.txt'.format(config['CHR'])),
        join(config["INPUT_DATA_DIR"], '1000GP_Phase3_chr{}.legend.gz'.format(config['CHR'])),
        join(config["INPUT_DATA_DIR"], '1000GP_Phase3_chr{}.hap.gz'.format(config['CHR'])),
        join(config["ARRAY_SNP_DIR"], "{}.bim".format(config["CHR"]))
    output:
        "out/0_raw/chr{}.snp".format(config['CHR']),
        "out/0_raw/chr{}.map".format(config['CHR']),
        "out/0_raw/chr{}.legend".format(config['CHR']),
        expand("out/0_raw/chr{}.{{pop}}.hap".format(config['CHR']), pop=config["POPS"])
    run:
        if not os.path.exists("out/0_raw"):
            os.makedirs("out/0_raw")
        assert (config['MAF_MODE'] in ['AND', 'OR']), 'MAF_MODE must be one of the [AND/OR]'
        extract_raw_data(raw_dir=config["INPUT_DATA_DIR"],
                         out_dir="out/0_raw", 
                         pops=config['POPS'],
                         chr_i=config['CHR'],
                         array_snp_file=join(config["ARRAY_SNP_DIR"], "{}.bim".format(config["CHR"])),
                         maf_threshold=config['MAF_THRESHOLD'],
                         maf_mode=config['MAF_MODE'])

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
        "out/1_hapgen/chr{}.EUR.phgeno".format(config['CHR']),
        "out/1_hapgen/chr{}.AFR.phgeno".format(config['CHR']),
        snp_file="out/0_raw/chr{}.snp".format(config['CHR']),
    output:
        "out/2_admix_geno/{prefix}/admix.phgeno"
    run:
        out_dir = output[0][0 : output[0].rfind('/')]
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        pop_prop, n_gen, n_sample = prefix2setting(wildcards.prefix)
        prefix_list = [f"out/1_hapgen/chr{config['CHR']}.{pop}" for pop in pop_prop.keys()]
        admix_prop = list(pop_prop.values()) 
	
        if int(n_sample) > config['HAPGEN_NUM_HAPLOS'] / 2:
            print("""
                  Warning: number of haplotypes in the admixture population should be smaller than 
                  half of ancestral population haplotypes. Please adjust numbers accordingly.
                  """)
        
        simulate_admixture(input_prefix_list=prefix_list, 
                           admix_prop=admix_prop,
                           num_gens=n_gen,
                           num_haplos=n_sample, 
                           snp_file=input.snp_file,
                           out_dir=out_dir, 
                           admix_simu_dir=config['ADMIX_SIMU_DIR'])

