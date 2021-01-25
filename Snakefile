from os.path import join
from utils import *
import os
import yaml

configfile: "config.yaml"

SIM_SETTINGS = [
    {"EUR": 0.5, "AFR": 0.5, "N_GEN": 10, "N_SAMPLE": 10000}
    ]

sim_code_list = [setting2prefix(*parse_sim_setting(sim)) for sim in SIM_SETTINGS]

rule all:
    input:
        "out/ukb_array/snp.txt",
        expand("out/ukb_array/{pop}.phgeno", pop=config["POPS"]),
        expand("out/ukb_array/{prefix}/admix.phgeno", prefix=sim_code_list)

rule extract_raw:
    input:
        join(config["KG_DATA_DIR"], 'genetic_map_chr{}_combined_b37.txt'.format(config['CHR'])),
        join(config["KG_DATA_DIR"], '1000GP_Phase3_chr{}.legend.gz'.format(config['CHR'])),
        join(config["KG_DATA_DIR"], '1000GP_Phase3_chr{}.hap.gz'.format(config['CHR'])),
        snp_list=join(config["DATA_DIR"], "{dataset}.snp_list")
    output:
        "out/{dataset}/snp.txt",
        "out/{dataset}/map.txt",
        "out/{dataset}/legend.txt",
        expand("out/{{dataset}}/{pop}.hap", pop=config["POPS"])
    run:
        out_dir = output[0].rsplit('/', 1)[0]
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        assert (config['MAF_MODE'] in ['AND', 'OR']), 'MAF_MODE must be one of the [AND/OR]'
        extract_raw_data(raw_dir=config["KG_DATA_DIR"],
                         out_dir=out_dir, 
                         pops=config['POPS'],
                         chr_i=config['CHR'],
                         snp_list=input.snp_list,
                         maf_threshold=config['MAF_THRESHOLD'],
                         maf_mode=config['MAF_MODE'])

rule extend_hapgen:
    input:
        map_file="{prefix}/map.txt",
        snp_file="{prefix}/snp.txt",
        legend_file="{prefix}/legend.txt",
        hap_file="{prefix}/{pop}.hap"
    output:
        "{prefix}/{pop}.phgeno"
    run:
        extend_hapgen(map_file=input.map_file,
                      legend_file=input.legend_file,
                      hap_file=input.hap_file,
                      num_haplos=config['HAPGEN_NUM_HAPLOS'],
                      hapgen2_path=config['HAPGEN2_PATH'],
                      out_prefix=output[0][0 : output[0].rfind('.')])
        

rule simulate_admixture:   
    input:
        "{dataset}/EUR.phgeno",
        "{dataset}/AFR.phgeno",
        snp_file="{dataset}/snp.txt"
    output:
        "{dataset}/{sim_code}/admix.phgeno"
    run:
        out_dir = output[0][0 : output[0].rfind('/')]
        pop_prop, n_gen, n_sample = prefix2setting(wildcards.sim_code)
        prefix_list = [f"{wildcards.dataset}/{pop}" for pop in pop_prop.keys()]
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

