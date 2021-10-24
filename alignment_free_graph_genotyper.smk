include:
    "alignment_free_graph_genotyper_prepare.smk"

def get_dataset_n_nodes(wildcards):
    return config["analysis_regions"][wildcards.dataset]["n_nodes"]

rule map:
    input:
        reads="data/{dataset}/{experiment}.fa",
        kmer_index_only_variants="data/{dataset}/kmer_index_only_variants.npz"
    output:
        node_counts="data/{dataset}/{experiment}.I{max_index_frequency}.node_counts.npy"
    params:
        n_nodes=get_dataset_n_nodes
    threads: config["n_threads"]
    benchmark:
        "data/{dataset}/benchmarks/mapI{max_index_frequency}_{experiment}.tsv"
    shell:
        "/usr/bin/time -v alignment_free_graph_genotyper count -i {input.kmer_index_only_variants} -n {output.node_counts} -t {config[n_threads]} -c 1000000 -r {input.reads} -M {params.n_nodes} --skip-chaining True -I {wildcards.max_index_frequency}"

def get_sample_name_from_experiment(wildcards):
    return wildcards.experiment.split("_")[0]


def get_read_coverage_from_experiment(wildcards):
    return wildcards.experiment.split("_")[-1].replace("x", "")

rule genotype:
    input:
        node_counts="data/{dataset}/{experiment}.I1000.node_counts.npy",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        variants="data/{dataset}/variants_no_overlaps_no_genotypes.vcf",
        model="data/{dataset}/combination_model.npz",
        genotype_frequencies="data/{dataset}/genotype_frequencies_{n_individuals}individuals.npz",
        most_similar_variant_lookup="data/{dataset}/most_similar_variant_lookup_{n_individuals}individuals.npz",
        transition_probs="data/{dataset}/transition_probs_{n_individuals}individuals.npy",
        tricky_variants="data/{dataset}/tricky_variants.npy",
    output:
        genotypes="data/{dataset}/usN{n_individuals,\d+}_{experiment,[a-z0-9_]+}.vcf.gz"
    threads:
        4
    benchmark:
        "data/{dataset}/benchmarks/usN{n_individuals}_{experiment}.tsv"
    params:
        sample_name=get_sample_name_from_experiment,
        read_coverage=get_read_coverage_from_experiment
    shell:
        "alignment_free_graph_genotyper genotype -c {input.node_counts} "
        "-g {input.variant_to_nodes} " 
        "-v {input.variants} " 
        "-A {input.model} " 
        "-G {input.genotype_frequencies} " 
        "-M {input.most_similar_variant_lookup} " 
        "-o {output.genotypes}.tmp " 
        "-C CombinationModelGenotyper " 
        "-t 10 -z 2000000 "
        "-a 2.0 -q -3 " 
        "-p {input.transition_probs} " 
        #"--average-coverage {params.read_coverage} "
        "--average-coverage 1 "
        "-x {input.tricky_variants} "
        "--sample-name-output {params.sample_name} "
        "&& bgzip -c {output.genotypes}.tmp > {output.genotypes} "


# hack to run genotyper with 2058 individuals if none are specified
rule genotype_wrapper:
    input:
        genotypes = "data/{dataset}/usN2058_{experiment}.vcf.gz"
    output:
        genotypes = "data/{dataset}/us_{experiment,[a-z0-9_]+}.vcf.gz"
    shell:
        "cp {input} {output}"


rule genotype_without_using_modelled_counts:
    input:
        node_counts="data/{dataset}/{experiment}.I1.node_counts.npy",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        variants="data/{dataset}/variants_no_overlaps_no_genotypes.vcf",
        genotype_frequencies="data/{dataset}/genotype_frequencies_naive.npz",
        tricky_variants="data/{dataset}/tricky_variants_nonunique_kmers.npy",
    output:
        genotypes="data/{dataset}/nomodel_us_{experiment}.vcf.gz"
    threads:
        4
    shell:
        "alignment_free_graph_genotyper genotype -c {input.node_counts} "
        "-g {input.variant_to_nodes} "
        "-v {input.variants} "
        "-G {input.genotype_frequencies} "
        "-o {output.genotypes}.tmp "
        "-C NumpyGenotyper "
        "-t 4 -z 500000 -a 1.0 -u True -z 6000 "
        "-x {input.tricky_variants}"
        "&& bgzip -c {output.genotypes}.tmp > {output.genotypes} "



