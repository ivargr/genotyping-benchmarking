
def get_dataset_n_nodes(wildcards):
    return config["analysis_regions"][wildcards.dataset]["n_nodes"]

include:
    "kage_prepare.smk"



rule map:
    input:
        reads="data/{dataset}/{experiment}.fa",
        kmer_index_only_variants="data/{dataset}/kmer_index_only_variants_with_revcomp.npz"
    output:
        node_counts="data/{dataset}/{experiment}.I{max_index_frequency}.node_counts.npy",
        benchmark_report="data/{dataset}/benchmarks/mapI{max_index_frequency}_{experiment}.tsv"
    params:
        n_nodes=get_dataset_n_nodes
    threads: config["n_threads"]
    shell:
        #"/usr/bin/time -v kage count -i {input.kmer_index_only_variants} -n {output.node_counts} -t {config[n_threads]} -c 1000000 -r {input.reads} -M {params.n_nodes} --skip-chaining True -I {wildcards.max_index_frequency}"
        "/usr/bin/time -v kmer_mapper map -i {input.kmer_index_only_variants} -o {output.node_counts} -t {config[n_threads]} -c 1250000 -f {input.reads} -I {wildcards.max_index_frequency} 2> {output.benchmark_report}"


def get_sample_name_from_experiment(wildcards):
    return wildcards.experiment.split("_")[0]


def get_read_coverage_from_experiment(wildcards):
    return wildcards.experiment.split("_")[-1].replace("x", "")

rule genotype:
    input:
        node_counts="data/{dataset}/{experiment}.I1000.node_counts.npy",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        #variants="data/{dataset}/variants_no_genotypes.vcf",
        variants="data/{dataset}/numpy_variants.npz",
        model="data/{dataset}/sampling_count_model_{n_individuals}all.npz",
        #genotype_frequencies="data/{dataset}/genotype_frequencies_{n_individuals}individuals.npz",
        #most_similar_variant_lookup="data/{dataset}/most_similar_variant_lookup_{n_individuals}individuals.npz",
        #transition_probs="data/{dataset}/transition_probs_{n_individuals}individuals.npy",
        helper_model="data/{dataset}/helper_model_{n_individuals}{subpopulation}.npy",
        helper_model_combo_matrix="data/{dataset}/helper_model_{n_individuals}{subpopulation}_combo_matrix.npy",
        tricky_variants="data/{dataset}/tricky_variants_{n_individuals}{subpopulation}.npy",
        #index_bundle="data/{dataset}/index_{n_individuals}{subpopulation}.npz"
    output:
        genotypes="data/{dataset}/usN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment,[A-Za-z0-9_]+}.vcf.gz",
        probs="data/{dataset}/usN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment,[A-Za-z0-9_]+}.vcf.gz.tmp.probs.npy",
        count_probs="data/{dataset}/usN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment,[A-Za-z0-9_]+}.vcf.gz.tmp.count_probs.npy",
        benchmark_report="data/{dataset}/benchmarks/usN{n_individuals}{subpopulation,[a-z]+}_{experiment}.tsv"
    threads:
        config["n_threads"]
    benchmark:
        "data/{dataset}/benchmarks/usN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment}-snakemake.tsv"
    conda: "envs/kage.yml"
    params:
        sample_name=get_sample_name_from_experiment,
        read_coverage=get_read_coverage_from_experiment
    shell:
        "/usr/bin/time -v kage genotype -c {input.node_counts} "
        #"-i {input.index_bundle} "
        "-g {input.variant_to_nodes} " 
        "-v {input.variants} " 
        "-A {input.model} " 
        "-f {input.helper_model} "
        "-F {input.helper_model_combo_matrix} "
        "-x {input.tricky_variants} "
        "-C CombinationModelGenotyper " 
        "-t {config[n_threads]} "
        "-q 0.5 " 
        "--average-coverage {params.read_coverage} "
        "-o {output.genotypes}.tmp " 
        "--sample-name-output {params.sample_name} 2> {output.benchmark_report} "
        "&& bgzip -c {output.genotypes}.tmp > {output.genotypes} "


rule genotype_without_helper_model:
    input:
        node_counts="data/{dataset}/{experiment}.I1000.node_counts.npy",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        #variants="data/{dataset}/variants_no_genotypes.vcf",
        variants="data/{dataset}/numpy_variants.npz",
        model="data/{dataset}/sampling_count_model_{n_individuals}all.npz",
        genotype_frequencies="data/{dataset}/genotype_frequencies_naive.npz",
        #most_similar_variant_lookup="data/{dataset}/most_similar_variant_lookup_{n_individuals}individuals.npz",
        #transition_probs="data/{dataset}/transition_probs_{n_individuals}individuals.npy",
        tricky_variants="data/{dataset}/tricky_variants_{n_individuals}{subpopulation}.npy",
        #index_bundle="data/{dataset}/index_{n_individuals}{subpopulation}.npz"
    output:
        genotypes="data/{dataset}/kageNoHelperModelN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment,[A-Za-z0-9_]+}.vcf.gz",
        probs="data/{dataset}/kageNoHelperModelN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment,[A-Za-z0-9_]+}.vcf.gz.tmp.probs.npy",
        count_probs="data/{dataset}/kageNoHelperModelN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment,[A-Za-z0-9_]+}.vcf.gz.tmp.count_probs.npy",
        #benchmark_report="data/{dataset}/benchmarks/kageNoHelperModelN{n_individuals}{subpopulation,[a-z]+}_{experiment}.tsv"
    threads:
        config["n_threads"]
    benchmark:
        "data/{dataset}/benchmarks/kageNoHelperModelN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment}-snakemake.tsv"
    conda: "envs/kage.yml"
    params:
        sample_name=get_sample_name_from_experiment,
        read_coverage=get_read_coverage_from_experiment
    shell:
        "/usr/bin/time -v kage genotype -c {input.node_counts} "
        #"-i {input.index_bundle} "
        "-g {input.variant_to_nodes} "
        "-v {input.variants} "
        "-A {input.model} "
        "-x {input.tricky_variants} "
        "-C CombinationModelGenotyper "
        "-G {input.genotype_frequencies} "
        "-t {config[n_threads]} "
        "-q 0.5 "
        "--average-coverage {params.read_coverage} "
        "-o {output.genotypes}.tmp "
        "--ignore-helper-variants True "
        "--sample-name-output {params.sample_name} "
        "&& bgzip -c {output.genotypes}.tmp > {output.genotypes} "


# hack to run genotyper with 2548 individuals if none are specified
rule genotype_wrapper:
    input:
        genotypes = "data/{dataset}/usN2548all_{experiment}.vcf.gz"
    output:
        genotypes = "data/{dataset}/us_{experiment,[a-z0-9_]+}.vcf.gz"
    shell:
        "cp {input} {output}"


rule kage_naive:
    input:
        node_counts="data/{dataset}/{experiment}.I1.node_counts.npy",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        variants="data/{dataset}/numpy_variants.npz",
        genotype_frequencies="data/{dataset}/genotype_frequencies_naive.npz",
        #tricky_variants="data/{dataset}/variants_with_nonunique_kmers.npy",
        tricky_variants="data/{dataset}/tricky_variants_2548all_nonunique_kmers.npy",  # NB! Using 2548 because we want variants thare are nonunique across whole population
        #tricky_variants="data/{dataset}/tricky_variants_2548all.npy",  # NB! Using 2548 because we want variants thare are nonunique across whole population
        #model="data/{dataset}/sampling_count_model_2548all.npz",
        helper_model="data/{dataset}/helper_model_{n_individuals}{subpopulation}.npy",
        helper_model_combo_matrix="data/{dataset}/helper_model_{n_individuals}{subpopulation}_combo_matrix.npy",
    output:
        genotypes="data/{dataset}/naivekageN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment}.vcf.gz"
    params:
        sample_name=get_sample_name_from_experiment,
        read_coverage=get_read_coverage_from_experiment
    threads:
        4
    shell:
        "kage genotype -c {input.node_counts} "
        "-g {input.variant_to_nodes} "
        "-v {input.variants} "
        "-G {input.genotype_frequencies} "
        "-o {output.genotypes}.tmp "
        #"-C NumpyGenotyper "
        #"-A {input.model} "
        "-C CombinationModelGenotyper "
        "--ignore-helper-variants True "
        "-t 1 -q 0.3 "
        "--average-coverage {params.read_coverage} "
        "--sample-name-output {params.sample_name} "
        "-f {input.helper_model} "
        "-F {input.helper_model_combo_matrix} "
        "-x {input.tricky_variants}"
        "&& bgzip -c {output.genotypes}.tmp > {output.genotypes} "



