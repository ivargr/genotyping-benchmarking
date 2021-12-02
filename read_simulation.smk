def get_simulation_tmp_datasets(wildcards):
    haplotypes = [0, 1]
    chromosomes = config["analysis_regions"][wildcards.dataset]["simulation_chromosomes"].split()
    files = []
    for chromosome in chromosomes:
        for haplotype in haplotypes:
            files.append("data/" + wildcards.dataset + "/" + wildcards.truth_dataset + "_raw_simulated_reads_chromosome" + chromosome + "_haplotype" + str(haplotype) + "_coverage" + wildcards.coverage + ".txt")

    return files


def get_real_data_reads_url(wildcards):
    return config["truth_datasets"][wildcards.truth_dataset]["reads_url"]


rule prepare_simulation:
    input:
        vcf="data/{dataset}/truth_{truth_dataset}.vcf.gz",
        reference="data/{dataset}/ref.fa"

    output:
        coordinate_map="data/{dataset}/{truth_dataset}_coordinate_map_chromosome{chromosome}_haplotype{haplotype}.npz",
        haplotype_reference="data/{dataset}/{truth_dataset}_chromosome{chromosome}_haplotype{haplotype}_reference.fasta",
        haplotype_reference_fasta="data/{dataset}/{truth_dataset}_chromosome{chromosome}_haplotype{haplotype}_reference.fasta.fai",
    conda: "envs/prepare_data.yml"
    shell:
        "graph_read_simulator prepare_simulation --chromosome {wildcards.chromosome} --haplotype {wildcards.haplotype} "
        "--vcf {input.vcf} --reference {input.reference} -o data/{wildcards.dataset}/{wildcards.truth_dataset}_ && "
        "samtools faidx {output.haplotype_reference} "


rule simulate_reads_for_chromosome_and_haplotype:
    input:
        coordinate_map="data/{dataset}/{truth_dataset}_coordinate_map_chromosome{chromosome}_haplotype{haplotype}.npz",
        haplotype_reference="data/{dataset}/{truth_dataset}_chromosome{chromosome}_haplotype{haplotype}_reference.fasta"

    output:
        "data/{dataset}/{truth_dataset}_raw_simulated_reads_chromosome{chromosome}_haplotype{haplotype}_coverage{coverage}.txt"
    conda: "envs/prepare_data.yml"
    shell:
        "graph_read_simulator simulate_reads -s 0.001 --deletion_prob 0.001 --insertion_prob 0.001 -D data/{wildcards.dataset}/{wildcards.truth_dataset}_ '{wildcards.chromosome} {wildcards.haplotype}' {wildcards.coverage} > {output}"

rule simulate_reads:
    input:
        get_simulation_tmp_datasets
    output:
        reads="data/{dataset}/{truth_dataset}_simulated_reads_{coverage,\d+}x.fa",
        read_positions="data/{dataset}/{truth_dataset}_simulated_reads_{coverage,\d+}x.readpositions"
    conda: "envs/prepare_data.yml"
    shell:
        "cat {input} | graph_read_simulator assign_ids {output.read_positions} {output.reads}"


rule get_real_raw_reads:
    output:
        reads="data/{dataset}/{truth_dataset}_real_reads_raw.fq.gz",
        #reads="data/{dataset}/{truth_dataset}_real_reads_{coverage,\d+}x.fq",
    params:
        url=get_real_data_reads_url
    conda: "envs/prepare_data.yml"
    shell:
        "wget -O - {params.url} | bedtools bamtofastq -i /dev/stdin -fq /dev/stdout | gzip > {output.reads}"


rule downsample_real_reads15x:
    input:
        rules.get_real_raw_reads.output
        #reads="data/{dataset}/{truth_dataset}_real_reads_raw.fq.gz",
    output:
        reads="data/{dataset}/{truth_dataset}_real_reads_15x.fq",
    conda: "envs/prepare_data.yml"
    shell:
        "zcat {input} | python scripts/downsample_fq.py 4 > {output.reads}"


rule downsample_real_reads30x:
    input:
        rules.get_real_raw_reads.output
        #reads="data/{dataset}/{truth_dataset}_real_reads_raw.fq.gz",
    output:
        reads="data/{dataset}/{truth_dataset}_real_reads_30x.fq",
    conda: "envs/prepare_data.yml"
    shell:
        "seqtk sample -2 -s{config[random_seed]} {input} 600000000 > {output.reads}"
        #"zcat {input} | python3 scripts/downsample_fq.py 4 > {output.reads}"


rule convert_real_reads_to_fa15x:
    input:
        rules.downsample_real_reads15x.output
        #"data/{dataset}/{truth_dataset}_real_reads_15x.fq",
    output:
        "data/{dataset}/{truth_dataset}_real_reads_15x.fa",
    conda: "envs/prepare_data.yml"
    shell:
        "cat {input} | seqtk seq -A > {output}"


rule convert_real_reads_to_fa30x:
    input:
        rules.downsample_real_reads30x.output
        #"data/{dataset}/{truth_dataset}_real_reads_15x.fq",
    output:
        "data/{dataset}/{truth_dataset}_real_reads_30x.fa",
    conda: "envs/prepare_data.yml"
    shell:
        "cat {input} | seqtk seq -A > {output}"

