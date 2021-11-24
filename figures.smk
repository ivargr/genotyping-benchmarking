
N_INDIVIDUALS_PANGENIE = [5, 15, 30, 50, 100]  #, 40, 50, 100, 200, 2058]
N_INDIVIDUALS=[5, 15, 30, 50, 100, 250, 500, 1000, 1750, 2548]  #, 40, 50, 100, 200, 2058]
#N_INDIVIDUALS=[5, 2058]  #, 40, 50, 100, 200, 2058]
WEB_FIGURE_DIR="/var/www/html/genotyping_figures/"

def figure2_file_names(wildcards):
    return ",".join(["data/dataset1/happy-hg002-usN" + str(i) + "_hg002_simulated_reads_15x-only-callable.summary.csv" for i in N_INDIVIDUALS] + \
         ["data/dataset1/happy-hg002-pangenieN" + str(i) + "_hg002_simulated_reads_15x-only-callable.summary.csv" for i in N_INDIVIDUALS_PANGENIE])

def figure2_names(wildcards):
    return ",".join(["us" for i in N_INDIVIDUALS] \
        + ["pangenie" for i in N_INDIVIDUALS_PANGENIE])

rule figure1:
    input:
        malva="data/dataset1/happy-hg002-malva_hg002_simulated_reads_15x.summary.csv",
        kage="data/dataset1/happy-hg002-usN2058_hg002_simulated_reads_15x.summary.csv",
        naivekage="data/dataset1/happy-hg002-naivekage_hg002_simulated_reads_15x.summary.csv"
    output:
        "figure1.html"
    shell:
        "genotyping_analysis plot_results_files -f {input.malva},{input.kage},{input.naivekage} -n malva,kage,naivekage -o {output}"


rule figure2:
    input:
        expand("data/dataset1/happy-hg002-usN{n_individuals}_hg002_simulated_reads_15x-only-callable.summary.csv", n_individuals=N_INDIVIDUALS),
        expand("data/dataset1/happy-hg002-pangenieN{n_individuals}_hg002_simulated_reads_15x-only-callable.summary.csv", n_individuals=N_INDIVIDUALS_PANGENIE),
        malva="data/dataset1/happy-hg002-malva_hg002_simulated_reads_15x-only-callable.summary.csv",
        naivekage="data/dataset1/happy-hg002-naivekage_hg002_simulated_reads_15x-only-callable.summary.csv"
    output:
        "figure2.html"
    params:
        file_names=figure2_file_names,
        names=figure2_names,
    shell:
        #"genotyping_analysis plot_results_files -f {input.malva},{input.us_no_model},{params.file_names} -n malva,nomodel,{params.names} -o {output}"
        "genotyping_analysis plot_results_files -f {params.file_names},{input.malva},{input.naivekage} -n {params.names},malva,naivekage -o {output} --type f1"


rule figure3:
    input:
        us = "data/dataset1/happy-hg002-us_hg002_simulated_reads_15x.summary.csv",
        graphtyper = "data/dataset1/happy-hg002-graphtyper_hg002_simulated_reads_15x.summary.csv",
        bayestyper = "data/dataset1/happy-hg002-bayestyper_hg002_simulated_reads_15x.summary.csv",
        gatk = "data/dataset1/happy-hg002-gatk_hg002_simulated_reads_15x.summary.csv",
        malva = "data/dataset1/happy-hg002-malva_hg002_simulated_reads_15x.summary.csv"
    output:
        "figure3.html"
    shell:
        "genotyping_analysis plot_results_files -f {input.us},{input.graphtyper},{input.bayestyper},{input.gatk} -n us,graphtyper,bayestyper,gatk -o {output}"


METHODS = ["us", "graphtyper", "bayestyper", "malva", "pangenie", "gatk"]
#METHODS = ["us", "graphtyper", "malva", "pangenie", "gatk"]
METHODS = ["us", "pangenie"]
METHODS_JOINED = ",".join(METHODS)


rule figure4:
    input:
        expand("data/dataset2/happy-hg002-{method}_hg002_real_reads_15x.summary.csv", method=METHODS)
    output:
        "figure4.html"
    run:
        file_names_joined = ",".join(input)
        shell("genotyping_analysis plot_results_files -f {file_names_joined} -n {METHODS_JOINED} -o {output}")


rule general_result_table:
    input:
        expand("data/{{dataset}}/happy-{{truth_dataset}}-{method}_{{experiment}}.summary.csv", method=METHODS)
    output:
        "table_{dataset,[a-z0-9_]+}-{experiment}-{truth_dataset,\w+}.html"
    shell:
        "python3 scripts/make_result_table.py {METHODS_JOINED} {wildcards.experiment} {wildcards.dataset} {wildcards.truth_dataset} > {output}"

rule simulated_data_result_table:
    input:
        "data/simulated_dataset2/happy-seed1-usN1000_seed1_simulated_reads_15x.summary.csv",
        "data/simulated_dataset2/happy-seed1-pangenie_seed1_simulated_reads_15x.summary.csv"
    output:
        "table_simulated_data.html"
    shell:
        "python3 scripts/make_result_table.py pangenie,usN1000 seed1_simulated_reads_15x simulated_dataset2 seed1 > {output}  && cat {output}"


rule table1:
    input:
        "table_dataset2-hg002_simulated_reads_15x-hg002.html"
        #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table1.html"
    shell:
        "cp {input} {output}"

rule table2:
    input:
        "table_dataset2-hg002_real_reads_15x-hg002.html"
    #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table2.html"
    shell:
        "cp {input} {output}"


rule table10:
    input:
        "table_dataset3-hg002_simulated_reads_15x-hg002.html"
    #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table10.html"
    shell:
        "cp {input} {output}"

rule table11:
    input:
        "table_dataset1-hg002_simulated_reads_15x-hg002.html"
    #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table11.html"
    shell:
        "cp {input} {output}"


rule table40:
    input:
        "table_dataset4-hg002_simulated_reads_15x-hg002.html"
    #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table40.html"
    shell:
        "cp {input} {output}"

rule table3:
    input:
        "table_dataset2-hg004_real_reads_15x-hg004.html"
    #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table3.html"
    shell:
        "cp {input} {output}"

rule table4:
    input:
        "table_dataset2-hg003_real_reads_15x-hg003.html"
    #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table4.html"
    shell:
        "cp {input} {output}"

rule move_figures:
    input:
        "{figure}.html"
    output:
        "/var/www/html/genotyping_figures/{figure}.html"
    shell:
        "cp {input} {output}"


rule trio_concordance_table:
    input:
        expand("data/dataset2/aj_pedigree_{method}/aj_pedigree_tab_delim_detailed_log.tsv", method=METHODS)




