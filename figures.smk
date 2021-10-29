
N_INDIVIDUALS_PANGENIE = [5, 15, 30, 50, 100]  #, 40, 50, 100, 200, 2058]
N_INDIVIDUALS=[5, 15, 30, 50, 100, 250, 500, 1000, 2058]  #, 40, 50, 100, 200, 2058]
N_INDIVIDUALS=[5, 2058]  #, 40, 50, 100, 200, 2058]
WEB_FIGURE_DIR="/var/www/html/genotyping_figures/"

def figure2_file_names(wildcards):
    return ",".join(["data/dataset1/happy-hg002-usN" + str(i) + "_hg002_simulated_reads_15x.summary.csv" for i in N_INDIVIDUALS] + \
         ["data/dataset1/happy-hg002-pangenie_hg002_simulated_reads_15x." + str(i) + "individuals.summary.csv" for i in N_INDIVIDUALS_PANGENIE])

def figure2_names(wildcards):
    return ",".join(["us" for i in N_INDIVIDUALS] \
        + ["pangenie" for i in N_INDIVIDUALS_PANGENIE])

rule figure1:
    input:
        malva="data/dataset1/happy-hg002-malva_hg002_simulated_reads_15x.summary.csv",
        us="data/dataset1/happy-hg002-usN2058_hg002_simulated_reads_15x.summary.csv",
        us_no_model="data/dataset1/happy-hg002-nomodel_us_hg002_simulated_reads_15x.summary.csv"
    output:
        "figure1.html"
    shell:
        "genotyping_analysis plot_results_files -f {input.malva},{input.us},{input.us_no_model} -n malva,us,nomodel -o {output}"


rule figure2:
    input:
        expand("data/dataset1/happy-hg002-usN{n_individuals}_hg002_simulated_reads_15x.summary.csv", n_individuals=N_INDIVIDUALS),
        expand("data/dataset1/happy-hg002-pangenie_hg002_simulated_reads_15x.{n_individuals}individuals.summary.csv", n_individuals=N_INDIVIDUALS_PANGENIE),
        malva="data/dataset1/happy-hg002-malva_hg002_simulated_reads_15x.summary.csv",
        #us_no_model="data/dataset1/happy-hg002-nomodel_us_hg002_simulated_reads_15x.summary.csv"
    output:
        "figure2.html"
    params:
        file_names=figure2_file_names,
        names=figure2_names,
    shell:
        #"genotyping_analysis plot_results_files -f {input.malva},{input.us_no_model},{params.file_names} -n malva,nomodel,{params.names} -o {output}"
        "genotyping_analysis plot_results_files -f {input.malva},{params.file_names} -n malva,{params.names} -o {output}"


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
METHODS = ["us", "graphtyper", "malva", "pangenie", "gatk"]
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


rule result_table:
    input:
        expand("data/{{dataset}}/happy-{{truth_dataset}}-{method}_{{experiment}}.summary.csv", method=METHODS)
    output:
        "table_{dataset,[a-z0-9]+}_{experiment}-{truth_dataset,\w+}.html"
    shell:
        "python3 scripts/make_result_table.py {METHODS_JOINED} {wildcards.experiment} {wildcards.dataset} {wildcards.truth_dataset} > {output}"


rule table1:
    input:
        "table_dataset2_hg002_simulated_reads_15x-hg002.html"
        #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table1.html"
    shell:
        "cp {input} {output}"

rule table2:
    input:
        "table_dataset2_hg002_real_reads_15x-hg002.html"
    #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table2.html"
    shell:
        "cp {input} {output}"

rule table3:
    input:
        "table_dataset2_hg004_real_reads_15x-hg004.html"
    #"table_dataset1_hg002_simulated_reads_15x-hg002.html"
    output:
        "table3.html"
    shell:
        "cp {input} {output}"

rule table4:
    input:
        "table_dataset2_hg003_real_reads_15x-hg003.html"
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




