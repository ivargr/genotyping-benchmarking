configfile: "config.yaml"

OUTPUT_PATH = ""

include:
    "analysis.smk"

include:
    "figures.smk"

include:
    "read_simulation.smk"

rule all:
    input:
        "data/simulated_dataset2/happy-seed1-usN50all_seed1_simulated_reads_15x.extended.csv"
        #OUTPUT_PATH + "table11.html"
