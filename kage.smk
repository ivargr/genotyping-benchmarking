include:
    "kage_prepare.smk"



def get_dataset_n_nodes(wildcards):
    return config["analysis_regions"][wildcards.dataset]["n_nodes"]

rule map:
    input:
        reads="data/{dataset}/{experiment}.fa",
        kmer_index_only_variants="data/{dataset}/kmer_index_only_variants.npz"
    output:
        node_counts="data/{dataset}/{experiment}.I{max_index_frequency}.node_counts.npy",
        benchmark_report="data/{dataset}/benchmarks/mapI{max_index_frequency}_{experiment}.tsv"
    params:
        n_nodes=get_dataset_n_nodes
    threads: config["n_threads"]
    #benchmark:
    #    "data/{dataset}/benchmarks/mapI{max_index_frequency}_{experiment}.tsv"
    shell:
        #"/usr/bin/time -v kage count -i {input.kmer_index_only_variants} -n {output.node_counts} -t {config[n_threads]} -c 1000000 -r {input.reads} -M {params.n_nodes} --skip-chaining True -I {wildcards.max_index_frequency}"
        "/usr/bin/time -v kmer_mapper map -i {input.kmer_index_only_variants} -o {output.node_counts} -t {config[n_threads]} -c 500000 -f {input.reads} -I {wildcards.max_index_frequency} 2> {output.benchmark_report}"


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
        model="data/{dataset}/combination_model.npz",
        #genotype_frequencies="data/{dataset}/genotype_frequencies_{n_individuals}individuals.npz",
        #most_similar_variant_lookup="data/{dataset}/most_similar_variant_lookup_{n_individuals}individuals.npz",
        #transition_probs="data/{dataset}/transition_probs_{n_individuals}individuals.npy",
        helper_model="data/{dataset}/helper_model_{n_individuals}individuals.npy",
        helper_model_combo_matrix="data/{dataset}/helper_model_{n_individuals}individuals_combo_matrix.npy",
        tricky_variants="data/{dataset}/tricky_variants.npy",
        index_bundle="data/{dataset}/index_{n_individuals}individuals.npz"
    output:
        genotypes="data/{dataset}/usN{n_individuals,\d+}_{experiment,[a-z0-9_]+}.vcf.gz",
        probs="data/{dataset}/usN{n_individuals,\d+}_{experiment,[a-z0-9_]+}.vcf.gz.tmp.probs.npy",
        count_probs="data/{dataset}/usN{n_individuals,\d+}_{experiment,[a-z0-9_]+}.vcf.gz.tmp.count_probs.npy",
        benchmark_report="data/{dataset}/benchmarks/usN{n_individuals}_{experiment}.tsv"
    threads:
        config["n_threads"]
    benchmark:
        "data/{dataset}/benchmarks/usN{n_individuals}_{experiment}-snakemake.tsv"
    conda: "envs/kage.yml"
    params:
        sample_name=get_sample_name_from_experiment,
        read_coverage=get_read_coverage_from_experiment
    shell:
        "/usr/bin/time -v kage genotype -c {input.node_counts} "
        "-i {input.index_bundle} "
        #"-g {input.variant_to_nodes} " 
        #"-v {input.variants} " 
        #"-A {input.model} " 
        #"-f {input.helper_model} "
        #"-F {input.helper_model_combo_matrix} "
        #"-x {input.tricky_variants} "
        "-C CombinationModelGenotyper " 
        "-t {config[n_threads]} "
        "-q 80 " 
        "--average-coverage {params.read_coverage} "
        "-o {output.genotypes}.tmp " 
        "--sample-name-output {params.sample_name} 2> {output.benchmark_report} "
        "&& bgzip -c {output.genotypes}.tmp > {output.genotypes} "



rule genotype_with_old_helper_model:
    input:
        node_counts="data/{dataset}/{experiment}.I1000.node_counts.npy",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        variants="data/{dataset}/variants_no_genotypes.vcf",
        model="data/{dataset}/combination_model.npz",
        genotype_frequencies="data/{dataset}/genotype_frequencies_{n_individuals}individuals.npz",
        most_similar_variant_lookup="data/{dataset}/most_similar_variant_lookup_{n_individuals}individuals.npz",
        transition_probs="data/{dataset}/transition_probs_{n_individuals}individuals.npy",
        #helper_model="data/{dataset}/helper_model_{n_individuals}individuals.npy",
        #helper_model_combo_matrix="data/{dataset}/helper_model_{n_individuals}individuals_combo_matrix.npy",
        tricky_variants="data/{dataset}/tricky_variants.npy",
    output:
        genotypes="data/{dataset}/kageoldN{n_individuals,\d+}_{experiment,[a-z0-9_]+}.vcf.gz",
        probs="data/{dataset}/kageoldN{n_individuals,\d+}_{experiment,[a-z0-9_]+}.vcf.gz.tmp.probs.npy"
    threads:
        4
    benchmark:
        "data/{dataset}/benchmarks/usN{n_individuals}_{experiment}.tsv"
    params:
        sample_name=get_sample_name_from_experiment,
        read_coverage=get_read_coverage_from_experiment
    conda: "envs/kage.yml"
    shell:
        "kage genotype -c {input.node_counts} "
        "-g {input.variant_to_nodes} "
        "-v {input.variants} "
        "-A {input.model} "
        "-G {input.genotype_frequencies} " 
        "-M {input.most_similar_variant_lookup} " 
        "-o {output.genotypes}.tmp "
        "-C NumpyGenotyper  "
        "-t 1 -z 2000000000000000 "
        "-q 80 "
        "-p {input.transition_probs} " 
        "--average-coverage 3.0 "
        "-x {input.tricky_variants} "
        "--sample-name-output {params.sample_name} "
        "&& bgzip -c {output.genotypes}.tmp > {output.genotypes} "



# hack to run genotyper with 2548 individuals if none are specified
rule genotype_wrapper:
    input:
        genotypes = "data/{dataset}/usN2548_{experiment}.vcf.gz"
    output:
        genotypes = "data/{dataset}/us_{experiment,[a-z0-9_]+}.vcf.gz"
    shell:
        "cp {input} {output}"


rule kage_naive:
    input:
        node_counts="data/{dataset}/{experiment}.I1.node_counts.npy",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        variants="data/{dataset}/numpy_variants.vcf",
        genotype_frequencies="data/{dataset}/genotype_frequencies_naive.npz",
        tricky_variants="data/{dataset}/tricky_variants_nonunique_kmers.npy",
    output:
        genotypes="data/{dataset}/naivekage_{experiment}.vcf.gz"
    threads:
        4
    shell:
        "kage genotype -c {input.node_counts} "
        "-g {input.variant_to_nodes} "
        "-v {input.variants} "
        "-G {input.genotype_frequencies} "
        "-o {output.genotypes}.tmp "
        "-C NumpyGenotyper "
        "-t 1  -a 1.0 -u True -z 60000000000 "
        "-x {input.tricky_variants}"
        "&& bgzip -c {output.genotypes}.tmp > {output.genotypes} "



