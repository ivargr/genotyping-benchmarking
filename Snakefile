configfile: "config.yaml"

WEB_PATH = "/var/www/html/genotyping_figures/"

include:
    "analysis.smk"

include:
    "figures.smk"

include:
    "read_simulation.smk"

include:
    "paragraph.smk"

include:
    "vg.smk"

include:
    "glimpse.smk"

include:
    "low_coverage_genotyping.smk"

include:
    "kage_with_mapped_reads.smk"

ruleorder:
    downsample_real_reads15x > convert_fa_to_fq

ruleorder:
    downsample_real_reads30x > convert_fa_to_fq


rule all:
    input:
        "data/dataset1/happy-hg002-kageNoHelperModelN250all_hg002_simulated_reads_15x.extended.csv",


