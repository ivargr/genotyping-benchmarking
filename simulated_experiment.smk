configfile: "config.yaml"

OUTPUT_PATH = ""

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
        OUTPUT_PATH + "table11.html"
