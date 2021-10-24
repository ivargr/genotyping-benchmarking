configfile: "config.yaml"

DATASETS = ["dataset1", "dataset2"]
WEB_PATH = "/var/www/html/genotyping_figures/"

include:
    "analysis.smk"

include:
    "figures.smk"

include:
    "read_simulation.smk"

ruleorder:
    downsample_real_reads > convert_fa_to_fq

rule all:
    input:
        WEB_PATH + "figure1.html"
        #WEB_PATH + "table4.html",
        #WEB_PATH + "table3.html",
        #WEB_PATH + "figure4.html",
        #WEB_PATH + "table2.html"
        #"/var/www/html/genotyping_figures/table1.html"
        #"/var/www/html/genotyping_figures/figure3.html"
        #"data/dataset1/hg002_real_reads_15x.fa"
        #"data/dataset1/gatk_hg002_simulated_reads_15x.vcf.gz"
        #"figure3.html"
        #"data/dataset1/bayestyper_hg002_simulated_reads_15x.vcf.gz"
        #"/var/www/html/genotyping_figures/figure3.html"
        #"data/dataset1/hg002_simulated_reads_1x.fa"
        #"data/dataset1/hg002_coordinate_map_chromosome1_haplotype0.npz"
        #"figure1.html"
        #"data/dataset1/happy-hg002-malva_simulated_reads15x.summary.csv"
        #"data/dataset1/malva_genotypes_simulated_reads15x.fa"
        #"data/dataset1/happy-hg002-genotypes_simulated_reads15x.summary.csv"

