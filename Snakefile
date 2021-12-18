configfile: "config.yaml"

# Change this to a path where you want figures to end up
WEB_PATH = "/var/www/html/genotyping_figures/"

include:
    "analysis.smk"

include:
    "figures.smk"

include:
    "read_simulation.smk"

ruleorder:
    downsample_real_reads15x > convert_fa_to_fq

ruleorder:
    downsample_real_reads30x > convert_fa_to_fq


rule all:
    input:
        WEB_PATH  + "figure1.html",  # Figure 3 in manuscript
        #WEB_PATH  + "figure2.html",  # Figure 4 in manuscript
        #WEB_PATH  + "table2.html",  # Table 2 in manuscript
        #WEB_PATH  + "table2.html",  # Table 2 in manuscript
        #WEB_PATH + "supplmentary_table.html",  # Supplementary table 1 in manuscript
        #WEB_PATH + "table20.html",  # Supplementary table 2 in manuscript
        #WEB_PATH + "supplementary_table2.html",  # Supplementary table 3 in manuscript

