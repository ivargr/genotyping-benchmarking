
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
        vcf = "data/{dataset}/variants_no_overlaps.vcf.gz",
        reference = "data/{dataset}/ref.fa",
    output:
        "data/{dataset}/obgraph_chr{chromosome}.npz"
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
    benchmark: "data/{dataset}/benchmarks/merge_chromosome_graphs.tsv"
    shell:
        "obgraph merge_graphs -o {output} -g {input} &&"
        "obgraph set_numeric_node_sequences -g {output}"


rule make_variant_to_nodes:
    input:
        graph="data/{dataset}/obgraph.npz",
        vcf="data/{dataset}/variants_no_overlaps_no_genotypes.vcf"
    output:
        "data/{dataset}/variant_to_nodes.npz"
    resources:
        mem_gb=60
    shell:
        "obgraph make_variant_to_nodes -g {input.graph} -v {input.vcf} -o {output}"


rule make_genotype_matrix:
    input:
        graph="data/{dataset}/obgraph.npz",
        vcf="data/{dataset}/variants_no_overlaps_{n_individuals}individuals.vcf.gz",
        vcf_no_genotypes="data/{dataset}/variants_no_overlaps_no_genotypes.vcf"
    output:
        genotype_matrix="data/{dataset}/genotype_matrix_{n_individuals}individuals.npy",
        most_similar_variant_lookup="data/{dataset}/most_similar_variant_lookup_{n_individuals}individuals.npz",
        transition_probs="data/{dataset}/transition_probs_{n_individuals}individuals.npy",
        genotype_frequencies="data/{dataset}/genotype_frequencies_{n_individuals}individuals.npz"
    threads: config["n_threads_data"]
    benchmark: "data/{dataset}/benchmarks/make_genotype_matrix_{n_individuals}individuals.tsv"
    resources:
        mem_gb=200
    shell:
        #"n_individuals=$(zcat {input.vcf} | head -n 1000 | grep -i '#chrom' | python3 -c 'import sys; print(len(list(sys.stdin)[0].split())-9)') || true && "
        "n_variants=$(grep -v '#' {input.vcf_no_genotypes} | wc -l) || true && "
        "obgraph make_genotype_matrix -t {config[n_threads_data]} -g {input.graph} -v {input.vcf} -o {output.genotype_matrix} -n {wildcards.n_individuals} -m $n_variants && "
        "obgraph analyse_genotype_matrix -t {config[n_threads_data]} -G {output.genotype_matrix} -o {output.most_similar_variant_lookup} &&"
        "obgraph make_transition_probabilities -t {config[n_threads_data]} -G {output.genotype_matrix} -m {output.most_similar_variant_lookup} -o {output.transition_probs} && "
        "obgraph get_genotype_frequencies -g {output.genotype_matrix} -o {output.genotype_frequencies} -t {config[n_threads_data]}"


rule make_naive_genotype_frequencies:
    input:
        vcf="data/{dataset}/variants_no_overlaps_no_genotypes.vcf"
    output:
        genotype_frequencies="data/{dataset}/genotype_frequencies_naive.npz"
    shell:
        "obgraph get_genotype_frequencies -v {input.vcf} -o {output} -t 40"


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
    benchmark:
        "data/{dataset}/benchmarks/sample_kmers_from_linear_reference.tsv"
    shell:
        "graph_kmer_index make -s 1 -k {config[k]} -G {params.genome_size} -o data/{wildcards.dataset}/linear_kmers -R {input} -n ref -t 20 && "
        "graph_kmer_index make_from_flat -o data/{wildcards.dataset}/linear_kmer_index -f {output.flat} -m 20000033"




rule make_variant_kmer_index:
    input:
        graph="data/{dataset}/obgraph.npz",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        vcf="data/{dataset}/variants_no_overlaps_no_genotypes.vcf",
        linear_kmer_index="data/{dataset}/linear_kmer_index.npz"
    output:
        variant_kmers="data/{dataset}/variant_kmers.npz",
        index="data/{dataset}/kmer_index_only_variants.npz",
        reverse_index="data/{dataset}/reverse_variant_kmers.npz"
    benchmark:
        "data/{dataset}/benchmarks/make_variant_kmer_index.tsv"
    shell:
        "graph_kmer_index make_unique_variant_kmers -g {input.graph} -V {input.variant_to_nodes} -k {config[k]} -o data/{wildcards.dataset}/variant_kmers -v {input.vcf} -t 30 -c 300000 --max-variant-nodes 5 -i {input.linear_kmer_index} && "
        "graph_kmer_index make_from_flat -o {output.index} -f data/{wildcards.dataset}/variant_kmers.npz -m 200000033 && "
        "graph_kmer_index make_reverse -f data/{wildcards.dataset}/variant_kmers -o {output.reverse_index}"


rule sample_kmers_for_kmer_model:
    input:
        graph="data/{dataset}/obgraph.npz",
        variant_kmers="data/{dataset}/variant_kmers.npz",
    params:
        genome_size=get_dataset_genome_size
    output:
        flat_kmers="data/{dataset}/flat_kmers_for_model.npz",
    benchmark: "data/{dataset}/benchmarks/sample_kmers_for_kmer_model.tsv"
    shell:
        "graph_kmer_index make -s 1 -k {config[k]} -G {params.genome_size} -o {output.flat_kmers} -g {input.graph} --max-frequency 100000 -v 4 -t 40 --only-save-one-node-per-kmer True --include-reverse-complement True --whitelist {input.variant_kmers}"


rule make_kmer_index:
    input:
        flat_kmers = "data/{dataset}/flat_kmers_for_model.npz",
    output:
        kmer_index = "data/{dataset}/kmer_index.npz"
    benchmark: "data/{dataset}/benchmarks/make_kmer_index.tsv"
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
    shell:
        "alignment_free_graph_genotyper model_using_kmer_index -i {input.kmer_index} -g {input.variant_to_nodes} -r {input.reverse_variant_kmers} -m {params.n_nodes} -o {output.haplotype_model} -k {config[k]} -t 8 &&"
        "alignment_free_graph_genotyper make_genotype_model -v {input.variant_to_nodes} -n {output.haplotype_model} -o {output.genotype_model}"
        
rule find_tricky_variants:
    input:
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        genotype_model= "data/{dataset}/genotype_model.npz",
        reverse_variant_kmers="data/{dataset}/reverse_variant_kmers.npz",
    output:
        tricky_variants="data/{dataset}/tricky_variants.npy",
        tricky_variants_nonunique_kmers="data/{dataset}/tricky_variants_nonunique_kmers.npy"
    shell:
        "alignment_free_graph_genotyper find_tricky_variants -v {input.variant_to_nodes} -m {input.genotype_model} -r {input.reverse_variant_kmers} -o {output.tricky_variants} -M 3 && "
        "alignment_free_graph_genotyper find_tricky_variants -v {input.variant_to_nodes} -m {input.genotype_model} -r {input.reverse_variant_kmers} -o {output.tricky_variants_nonunique_kmers} -u True"

