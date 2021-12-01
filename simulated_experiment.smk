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
        OUTPUT_PATH + "table11.html"
