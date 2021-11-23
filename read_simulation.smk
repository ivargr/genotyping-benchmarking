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

    shell:
        "graph_read_simulator simulate_reads -s 0.001 --deletion_prob 0.001 --insertion_prob 0.001 -D data/{wildcards.dataset}/{wildcards.truth_dataset}_ '{wildcards.chromosome} {wildcards.haplotype}' {wildcards.coverage} > {output}"

rule simulate_reads:
    input:
        get_simulation_tmp_datasets
    output:
        reads="data/{dataset}/{truth_dataset}_simulated_reads_{coverage,\d+}x.fa",
        read_positions="data/{dataset}/{truth_dataset}_simulated_reads_{coverage,\d+}x.readpositions"

    shell:
        "cat {input} | graph_read_simulator assign_ids {output.read_positions} {output.reads}"


rule get_real_raw_reads:
    output:
        reads="data/{dataset}/{truth_dataset}_real_reads_raw.fq.gz",
        #reads="data/{dataset}/{truth_dataset}_real_reads_{coverage,\d+}x.fq",
    params:
        url=get_real_data_reads_url
    shell:
        "wget -O - {params.url} | bedtools bamtofastq -i /dev/stdin -fq /dev/stdout | gzip > {output.reads}"


rule downsample_real_reads:
    input:
        rules.get_real_raw_reads.output
        #reads="data/{dataset}/{truth_dataset}_real_reads_raw.fq.gz",
    output:
        reads="data/{dataset}/{truth_dataset}_real_reads_15x.fq",
    shell:
        "zcat {input} | python3 scripts/downsample_fq.py 4 > {output.reads}"


rule convert_real_reads_to_fa:
    input:
        rules.downsample_real_reads.output
        #"data/{dataset}/{truth_dataset}_real_reads_15x.fq",
    output:
        "data/{dataset}/{truth_dataset}_real_reads_15x.fa",
    shell:
        "cat {input} | seqtk seq -A > {output}"

