
N_INDIVIDUALS_PANGENIE = [5, 15, 30, 50, 100]  #, 40, 50, 100, 200, 2058]
N_INDIVIDUALS=[5, 15, 30, 50, 100, 250, 500, 1000, 1750, 2548]  #, 40, 50, 100, 200, 2058]
#N_INDIVIDUALS=[5, 15, 30, 50, 100, 250, 500, 550, 600, 650, 700, 800, 1000, 2548]  #, 40, 50, 100, 200, 2058]

METHODS = ["us", "graphtyper", "bayestyper", "malva", "pangenie", "gatk"]
#METHODS = ["us", "graphtyper", "malva", "pangenie", "gatk"]
METHODS = ["us", "pangenie", "bayestyper", "malva", "graphtyper", "gatk"]
#METHODS = ["usN150all", "paragraph"]

METHODS_JOINED = ",".join(METHODS)

def figure2_file_names(wildcards):
    return ",".join(["data/dataset1/happy-hg002-usN" + str(i) + "all_hg002_simulated_reads_15x.extended.csv" for i in N_INDIVIDUALS] + \
         ["data/dataset1/happy-hg002-pangenieN" + str(i) + "_hg002_simulated_reads_15x.extended.csv" for i in N_INDIVIDUALS_PANGENIE])

def figure2_names(wildcards):
    return ",".join(["us" for i in N_INDIVIDUALS] \
        + ["pangenie" for i in N_INDIVIDUALS_PANGENIE])

rule figure1:
    input:
        malva="data/dataset1/happy-hg002-malva_hg002_simulated_reads_15x.extended.csv",
        kage="data/dataset1/happy-hg002-kageNoHelperModelN250all_hg002_simulated_reads_15x.extended.csv",
        naivekage="data/dataset1/happy-hg002-naivekageN250all_hg002_simulated_reads_15x.extended.csv"
    output:
        "figure1.html"
    conda: "envs/analysis.yml"
    shell:
        "python scripts/make_scatter_plot.py plot_results_files -f {input.malva},{input.kage},{input.naivekage} -n malva,kageNoHelperModel,naivekage -o {output}"


rule figure2:
    input:
        expand("data/dataset1/happy-hg002-usN{n_individuals}all_hg002_simulated_reads_15x.extended.csv", n_individuals=N_INDIVIDUALS),
        expand("data/dataset1/happy-hg002-pangenieN{n_individuals}_hg002_simulated_reads_15x.extended.csv", n_individuals=N_INDIVIDUALS_PANGENIE),
        malva="data/dataset1/happy-hg002-malva_hg002_simulated_reads_15x.extended.csv",
        naivekage="data/dataset1/happy-hg002-naivekage_hg002_simulated_reads_15x.extended.csv"
    output:
        "figure2.html"
    params:
        file_names=figure2_file_names,
        names=figure2_names,
    conda: "envs/analysis.yml"
    shell:
        #"python3 scripts/make_scatter_plot.py plot_results_files -f {input.malva},{input.us_no_model},{params.file_names} -n malva,nomodel,{params.names} -o {output}"
        "python scripts/make_scatter_plot.py plot_results_files -f {params.file_names},{input.malva},{input.naivekage} -n {params.names},malva,naivekage -o {output} --type f1"


rule figure3:
    input:
        us = "data/dataset1/happy-hg002-us_hg002_simulated_reads_15x.extended.csv",
        graphtyper = "data/dataset1/happy-hg002-graphtyper_hg002_simulated_reads_15x.extended.csv",
        bayestyper = "data/dataset1/happy-hg002-bayestyper_hg002_simulated_reads_15x.extended.csv",
        gatk = "data/dataset1/happy-hg002-gatk_hg002_simulated_reads_15x.extended.csv",
        malva = "data/dataset1/happy-hg002-malva_hg002_simulated_reads_15x.extended.csv"
    output:
        "figure3.html"
    conda: "envs/analysis.yml"
    shell:
        "python scripts/make_scatter_plot.py plot_results_files -f {input.us},{input.graphtyper},{input.bayestyper},{input.gatk} -n us,graphtyper,bayestyper,gatk -o {output}"



rule figure4:
    input:
        expand("data/dataset2/happy-hg002-{method}_hg002_real_reads_15x.extended.csv", method=METHODS)
    output:
        "figure4.html"
    run:
        file_names_joined = ",".join(input)
        shell("python3 scripts/make_scatter_plot.py plot_results_files -f {file_names_joined} -n {METHODS_JOINED} -o {output}")


rule general_result_table:
    input:
        expand("data/{{dataset}}/happy-{{truth_dataset}}-{method}_{{experiment}}.extended.csv", method=METHODS)
    output:
        table="table_{dataset,[a-z0-9_]+}-{experiment}-{truth_dataset,\w+}.html",
        f1_figure="f1_{dataset,[a-z0-9_]+}-{experiment}-{truth_dataset,\w+}.html",
        mem_figure="memory_{dataset,[a-z0-9_]+}-{experiment}-{truth_dataset,\w+}.html",
        runtime_figure="runtime_{dataset,[a-z0-9_]+}-{experiment}-{truth_dataset,\w+}.html",
    conda: "envs/analysis.yml"
    shell:
        "python scripts/make_result_table.py {METHODS_JOINED} {wildcards.experiment} {wildcards.dataset} {wildcards.truth_dataset} _{wildcards.dataset}-{wildcards.experiment}-{wildcards.truth_dataset}.html > {output.table}"

rule simulated_data_result_table:
    input:
        "data/simulated_dataset2/happy-seed1-usN1000all_seed1_simulated_reads_15x.extended.csv",
        "data/simulated_dataset2/happy-seed1-pangenie_seed1_simulated_reads_15x.extended.csv"
    output:
        "table_simulated_data.html"
    conda: "envs/analysis.yml"
    shell:
        "python scripts/make_result_table.py pangenie,usN1000 seed1_simulated_reads_15x simulated_dataset2 seed1 > {output}  && cat {output}"


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


rule table20:
    input:
        "table_dataset2-hg006_real_reads_15x-hg006.html"
    #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table20.html"
    shell:
        "cp {input} {output}"


rule supplementary_table:
    input:
        "table_dataset2-hg002_real_reads_30x-hg002.html"
    output:
        "supplementary_table.html"
    shell:
        "cp {input} {output}"

rule supplementary_table2:
    input:
        "table_dataset2-hg006_real_reads_30x-hg006.html"
    output:
        "supplementary_table2.html"
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
        "table_simulated_dataset0-seed1_simulated_reads_10x-seed1.html"
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
    
rule svtable1:
    input:
        "table_svdataset2-NA12889_simulated_reads_15x-NA12889.html"
    output:
        "svtable1.html"
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


rule f1_figure:
    input:
        expand("data/{{dataset}}/happy-{{truth_dataset}}-{method}_{{experiment}}.extended.csv",method=METHODS)
    output:
        "f1_figure_{dataset,[a-z0-9_]+}-{experiment}-{truth_dataset,\w+}.html"
    conda: "envs/analysis.yml"
    shell:
        "python scripts/make_f1_figure.py {METHODS_JOINED} {wildcards.experiment} {wildcards.dataset} {wildcards.truth_dataset} {output}"


rule f1_figure_whole_genome:
    input:
        "f1_figure_dataset2-hg002_real_reads_15x-hg002.html"
    output:
        "f1_figure_whole_genome.html"
    shell:
        "cp {input} {output}"

