
include:
    "prepare_data.smk"


def get_dataset_chromosomes(wildcards):
    return config["analysis_regions"][wildcards.dataset]["chromosomes"]

def get_dataset_n_nodes(wildcards):
    return config["analysis_regions"][wildcards.dataset]["n_nodes"]

def get_dataset_chromosomes_as_list(wildcards):
    return config["analysis_regions"][wildcards.dataset]["chromosomes"].split()

def get_chromosome_graph_names(wildcards):
    return ["data/" + wildcards.dataset + "/obgraph_chr%s.npz"  % chromosome for chromosome in config["analysis_regions"][wildcards.dataset]["chromosomes"].split()]

def get_dataset_genome_size(wildcards):
    return config["analysis_regions"][wildcards.dataset]["genome_size"]


rule make_chromosome_graph:
    input:
        vcf = "data/{dataset}/variants.vcf.gz",
        reference = "data/{dataset}/ref.fa",
    output:
        "data/{dataset}/obgraph_chr{chromosome}.npz"
    conda: "envs/kage.yml"
    shell:
        "obgraph make -r {input.reference} -c {wildcards.chromosome} -v {input.vcf} -o {output} && "
        "obgraph add_allele_frequencies -c {wildcards.chromosome} -v {input.vcf} -g {output}"


rule merge_chromosome_graphs:
    input:
        get_chromosome_graph_names
    output:
        "data/{dataset}/obgraph.npz"
    params:
        chromosomes=get_dataset_chromosomes
    resources:
        mem_gb=80
    conda: "envs/kage.yml"
    benchmark: "data/{dataset}/benchmarks/merge_chromosome_graphs.tsv"
    shell:
        "obgraph merge_graphs -o {output} -g {input} &&"
        "obgraph set_numeric_node_sequences -g {output}"


rule make_variant_to_nodes:
    input:
        graph="data/{dataset}/obgraph.npz",
        vcf="data/{dataset}/variants_no_genotypes.vcf"
    output:
        "data/{dataset}/variant_to_nodes.npz"
    resources:
        mem_gb=60
    conda: "envs/kage.yml"
    shell:
        "obgraph make_variant_to_nodes -g {input.graph} -v {input.vcf} -o {output}"


rule make_genotype_matrix:
    input:
        graph="data/{dataset}/obgraph.npz",
        vcf="data/{dataset}/variants_{n_individuals}individuals.vcf.gz",
        vcf_no_genotypes="data/{dataset}/variants_no_genotypes.vcf"
    output:
        genotype_matrix="data/{dataset}/genotype_matrix_{n_individuals,\d+}individuals.npy",
    threads: config["n_threads_data_quarter"]
    benchmark: "data/{dataset}/benchmarks/make_genotype_matrix_{n_individuals}individuals.tsv"
    resources:
        mem_gb=230
    conda: "envs/kage.yml"
    shell:
        #"n_individuals=$(zcat {input.vcf} | head -n 1000 | grep -i '#chrom' | python3 -c 'import sys; print(len(list(sys.stdin)[0].split())-9)') || true && "
        "n_variants=$(grep -v '#' {input.vcf_no_genotypes} | wc -l) || true && "
        "obgraph make_genotype_matrix -g {input.graph} -v {input.vcf} -o {output.genotype_matrix} -n {wildcards.n_individuals} -m $n_variants -t {config[n_threads_data]}"


"""
rule convert_genotype_matrix:
    input:
        genotype_matrix = "data/{dataset}/genotype_matrix_{n_individuals}individuals.npy",
    output:
        genotype_matrix = "data/{dataset}/genotype_matrix_converted_{n_individuals,\d+}individuals.npy",
    resources:
        mem_gb=200
    threads:
        config["n_threads_data"]
    shell:
        "obgraph convert_genotype_matrix -g {input.genotype_matrix} -o {output.genotype_matrix} -t 8"
"""

rule make_transition_probabilities:
    input:
        genotype_matrix="data/{dataset}/genotype_matrix_{n_individuals}individuals.npy",
        most_similar_variant_lookup="data/{dataset}/most_similar_variant_lookup_{n_individuals}individuals.npz",
    output:
        transition_probs = "data/{dataset}/transition_probs_{n_individuals}individuals.npy",
    threads: config["n_threads_data_quarter"]
    benchmark: "data/{dataset}/benchmarks/make_transition_probabilities_{n_individuals}individuals.tsv"
    resources: mem_gb=100
    conda: "envs/kage.yml"
    shell:
        "obgraph make_transition_probabilities -t {config[n_threads_data_quarter]} -G {input.genotype_matrix} -m {input.most_similar_variant_lookup} -o {output.transition_probs}"


rule make_most_similar_variant_lookup:
    input:
        genotype_matrix="data/{dataset}/genotype_matrix_{n_individuals}individuals.npy",
    output:
        most_similar_variant_lookup="data/{dataset}/most_similar_variant_lookup_{n_individuals}individuals.npz",
    threads: config["n_threads_data_quarter"]
    benchmark: "data/{dataset}/benchmarks/make_most_similar_variant_lookup_{n_individuals}individuals.tsv"
    resources: mem_gb=100
    conda: "envs/kage.yml"
    shell:
        "obgraph analyse_genotype_matrix -t {config[n_threads_data_quarter]} -G {input.genotype_matrix} -o {output.most_similar_variant_lookup}"


rule get_genotype_frequencies:
    input:
        genotype_matrix = "data/{dataset}/genotype_matrix_{n_individuals}individuals.npy",
    output:
        genotype_frequencies="data/{dataset}/genotype_frequencies_{n_individuals}individuals.npz"
    threads: config["n_threads_data"]
    benchmark: "data/{dataset}/benchmarks/get_genotype_frequencies_{n_individuals}individuals.tsv"
    resources:
        mem_gb=100
    conda: "envs/kage.yml"
    shell:
        "obgraph get_genotype_frequencies -g {input.genotype_matrix} -o {output.genotype_frequencies} -t {config[n_threads_data]}"


rule make_naive_genotype_frequencies:
    input:
        vcf="data/{dataset}/variants_no_genotypes.vcf"
    output:
        genotype_frequencies="data/{dataset}/genotype_frequencies_naive.npz"
    threads: config["n_threads_data"]
    benchmark: "data/{dataset}/benchmarks/get_naive_genotype_frequencies.tsv"
    resources:
        mem_gb = 60
    conda: "envs/kage.yml"
    shell:
        "obgraph get_genotype_frequencies -v {input.vcf} -o {output} -t {config[n_threads_data]}"


rule sample_kmers_from_linear_reference:
    input:
        "data/{dataset}/ref_flat.fa"
    output:
        flat="data/{dataset}/linear_kmers.npz",
        index="data/{dataset}/linear_kmer_index.npz"
    params:
        genome_size=get_dataset_genome_size
    resources:
        mem_gb=240
    threads: config["n_threads_data"]
    benchmark:
        "data/{dataset}/benchmarks/sample_kmers_from_linear_reference.tsv"
    conda: "envs/kage.yml"
    shell:
        "graph_kmer_index make -s 1 -k {config[k]} -G {params.genome_size} -o data/{wildcards.dataset}/linear_kmers -R {input} -n ref -t {config[n_threads_data]} && "
        "graph_kmer_index make_from_flat -o data/{wildcards.dataset}/linear_kmer_index -f {output.flat} -m 20000033"


rule make_helper_model:
    input:
        genotype_matrix="data/{dataset}/genotype_matrix_{n_individuals}individuals.npy",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        node_count_model="data/{dataset}/combination_model.npz",
    output:
        helper_model="data/{dataset}/helper_model_{n_individuals}individuals.npy",
        helper_model_combo_matrix="data/{dataset}/helper_model_{n_individuals}individuals_combo_matrix.npy"
    threads: config["n_threads_data"]
    resources: mem_gb=300
    conda: "envs/kage.yml"
    shell:
        "kage create_helper_model -o data/{wildcards.dataset}/helper_model_{wildcards.n_individuals}individuals "
        "-g {input.genotype_matrix} -w 100 -v {input.variant_to_nodes} -n {input.node_count_model} -t {config[n_threads_data]} "
        " "

rule make_variant_kmer_index:
    input:
        graph="data/{dataset}/obgraph.npz",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        vcf="data/{dataset}/variants_no_genotypes.vcf",
        linear_kmer_index="data/{dataset}/linear_kmer_index.npz"
    output:
        variant_kmers="data/{dataset}/variant_kmers.npz",
        index="data/{dataset}/kmer_index_only_variants.npz",
        reverse_index="data/{dataset}/reverse_variant_kmers.npz"
    benchmark:
        "data/{dataset}/benchmarks/make_variant_kmer_index.tsv"
    resources: mem_gb=100
    threads: config["n_threads_data"]
    conda: "envs/kage.yml"
    shell:
        "graph_kmer_index make_unique_variant_kmers -g {input.graph} -V {input.variant_to_nodes} -k {config[k]} -o data/{wildcards.dataset}/variant_kmers -v {input.vcf} -t {config[n_threads_data]} -c 10000 --max-variant-nodes 5 -i {input.linear_kmer_index} && "
        "graph_kmer_index make_from_flat -o {output.index} -f data/{wildcards.dataset}/variant_kmers.npz -m 200000033 && "
        "graph_kmer_index make_reverse -f data/{wildcards.dataset}/variant_kmers -o {output.reverse_index}"


rule make_numpy_variants:
    input:
        "data/{dataset}/variants_no_genotypes.vcf"
    output:
        "data/{dataset}/numpy_variants.npz"
    conda: "envs/kage.yml"
    shell:
        "obgraph make_numpy_variants -v {input} -o {output}"


rule sample_kmers_for_kmer_model:
    input:
        graph="data/{dataset}/obgraph.npz",
        variant_kmers="data/{dataset}/variant_kmers.npz",
    params:
        genome_size=get_dataset_genome_size
    output:
        flat_kmers="data/{dataset}/flat_kmers_for_model.npz",
    benchmark: "data/{dataset}/benchmarks/sample_kmers_for_kmer_model.tsv"
    resources: mem_gb=100
    threads: config["n_threads_data"]
    conda: "envs/kage.yml"
    shell:
        "graph_kmer_index make -s 1 -k {config[k]} -G {params.genome_size} -o {output.flat_kmers} -g {input.graph} --max-frequency 100000 -v 4 -t {config[n_threads_data]} --only-save-one-node-per-kmer True --include-reverse-complement True --whitelist {input.variant_kmers}"


rule make_kmer_index:
    input:
        flat_kmers = "data/{dataset}/flat_kmers_for_model.npz",
    output:
        kmer_index = "data/{dataset}/kmer_index.npz"
    benchmark: "data/{dataset}/benchmarks/make_kmer_index.tsv"
    resources:
        mem_gb=20
    conda: "envs/kage.yml"
    shell:
        "graph_kmer_index make_from_flat -o {output.kmer_index} -f {input.flat_kmers} -m 20000033 "


rule make_kmer_count_model:
    input:
        kmer_index="data/{dataset}/kmer_index.npz",
        reverse_variant_kmers="data/{dataset}/reverse_variant_kmers.npz",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz"
    output:
        haplotype_model="data/{dataset}/model.npz",
        genotype_model= "data/{dataset}/genotype_model.npz",
    benchmark: "data/{dataset}/benchmarks/make_kmer_count_model.tsv"
    params:
        n_nodes=get_dataset_n_nodes
    threads: config["n_threads_data_quarter"]
    resources:
        mem_gb=120
    conda: "envs/kage.yml"
    shell:
        "kage model -i {input.kmer_index} -g {input.variant_to_nodes} -r {input.reverse_variant_kmers} -m {params.n_nodes} -o {output.haplotype_model} -k {config[k]} -t {config[n_threads_data_quarter]} -V '' &&"
        "kage make_genotype_model -v {input.variant_to_nodes} -n {output.haplotype_model} -o {output.genotype_model}"

rule make_combination_model:
    input:
        kmer_index="data/{dataset}/kmer_index.npz",
        reverse_variant_kmers="data/{dataset}/reverse_variant_kmers.npz",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz"
    output:
        model="data/{dataset}/combination_model.npz",
    benchmark: "data/{dataset}/benchmarks/make_combination_model.tsv"
    params:
        n_nodes=get_dataset_n_nodes
    threads: config["n_threads_data_quarter"]
    resources:
        mem_gb=120
    conda: "envs/kage.yml"
    shell:
        "kage model -i {input.kmer_index} -g {input.variant_to_nodes} -r "
        "{input.reverse_variant_kmers} -m {params.n_nodes} -o {output.model} "
        "-k {config[k]} -t {config[n_threads_data_quarter]} -V v3 "


rule find_tricky_variants:
    input:
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        model= "data/{dataset}/combination_model.npz",
        reverse_variant_kmers="data/{dataset}/reverse_variant_kmers.npz",
    output:
        tricky_variants="data/{dataset}/tricky_variants.npy",
        tricky_variants_nonunique_kmers="data/{dataset}/tricky_variants_nonunique_kmers.npy"
    conda: "envs/kage.yml"
    shell:
        "kage find_tricky_variants -v {input.variant_to_nodes} -m {input.model} -r {input.reverse_variant_kmers} -o {output.tricky_variants} -M 1000 && "
        "kage find_tricky_variants -v {input.variant_to_nodes} -m {input.model} -r {input.reverse_variant_kmers} -o {output.tricky_variants_nonunique_kmers} -u True"

