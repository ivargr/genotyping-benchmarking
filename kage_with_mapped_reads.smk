
def get_dataset_n_nodes(wildcards):
    return config["analysis_regions"][wildcards.dataset]["n_nodes"]

rule vg_giraffe_map:
    input:
        xg="data/{dataset}/graph.xg",
        dist = "data/{dataset}/graph.dist",
        min = "data/{dataset}/graph.min",
        gbwt = "data/{dataset}/graph_downsampled.gbwt",
        gg = "data/{dataset}/graph_downsampled.gg",
        fq = "data/{dataset}/{reads}.fq"
    output:
        "data/{dataset}/{reads}.gaf"
    shell:
        "vg giraffe --hit-cap 100 --hard-hit-cap 5000 --max-multimaps 10000 -t {config[n_threads]} -p -m {input.min} -d {input.dist} -x {input.xg} -g {input.gg} -H {input.gbwt} -o gaf -f {input.fq} > {output}"

def get_sample_name_from_experiment(wildcards):
    return wildcards.experiment.split("_")[0]


def get_read_coverage_from_experiment(wildcards):
    return wildcards.experiment.split("_")[-1].replace("x", "")


rule get_node_counts_from_giraffe_mapping:
    input:
        gaf="data/{dataset}/{reads}.gaf",
        edge_mapping="data/{dataset}/obgraph.npz.edge_mapping",
        graph="data/{dataset}/obgraph.npz"
    output:
        "data/{dataset}/{reads}.node_counts_giraffe.npy"
    params:
        n_nodes = get_dataset_n_nodes
    shell:
        "kage node_counts_from_gaf -g {input.gaf} --min-mapq 30 --min-score 135 -o {output} -m {input.edge_mapping} -G {input.graph}"


rule simulate_reads_to_make_model_for_mapped_reads:
    input:
        ref="data/{dataset}/ref.fa"
    output:
        sampled_reads="data/{dataset}/sampled_reads.fa",
        positions="data/{dataset}/sampled_reads.txt"
    shell:
        "graph_read_simulator sample_from_linear_ref -f {input.ref} -s 10 | graph_read_simulator assign_ids {output.positions} {output.sampled_reads}"


rule store_sampled_reads_as_numpy_alignments:
    input:
        "data/{dataset}/sampled_reads.txt"
    output:
        "data/{dataset}/sampled_reads.npz"
    shell:
        "n_reads=$(cat {input} | wc -l || true) && "
        "cat {input} | numpy_alignments store truth {output} $n_reads"


rule giraffe_map_sampled_model_reads:
    input:
        xg = "data/{dataset}/graph.xg",
        dist = "data/{dataset}/graph.dist",
        min = "data/{dataset}/graph.min",
        gbwt = "data/{dataset}/graph_downsampled.gbwt",
        gg = "data/{dataset}/graph_downsampled.gg",
        fq = "data/{dataset}/sampled_reads.fa"
    output:
        "data/{dataset}/sampled.gaf"
    shell:
        "vg giraffe --hit-cap 100 --hard-hit-cap 5000 --max-multimaps 10000 -t {config[n_threads_data]} -p -m {input.min} -d {input.dist} -x {input.xg} -g "
        "{input.gg} -H {input.gbwt} -o gaf -f {input.fq} > {output}"
    
    
rule make_kage_model_for_mapped_reads:
    input:
        graph="data/{dataset}/obgraph.npz",
        sampled_read_positions="data/{dataset}/sampled_reads.npz",
        edge_mapping="data/{dataset}/obgraph.npz.edge_mapping",
        gaf="data/{dataset}/sampled.gaf"
    output:
        "data/{dataset}/model_for_read_mapping_with_sampled_reads.npz"
    shell:
        "kage model_for_read_mapping_using_sampled_reads -g {input.graph} -r {input.gaf} "
        "-p {input.sampled_read_positions} -o {output} -m 135 -e {input.edge_mapping} -c 15"
    

rule make_naive_kage_model_for_mapped_reads:
    input:
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz"
    output:
        "data/{dataset}/model_for_read_mapping.npz"
    shell:
        "kage model_for_read_mapping -g {input.variant_to_nodes} -o {output} "


rule kage_with_mapped_reads:
    input:
        node_counts="data/{dataset}/{experiment}.node_counts_giraffe.npy",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        variants="data/{dataset}/numpy_variants.npz",
        #model="data/{dataset}/model_for_read_mapping_with_sampled_reads.npz",
        model="data/{dataset}/model_for_read_mapping.npz",
        helper_model="data/{dataset}/helper_model_{n_individuals}{subpopulation}.npy",
        helper_model_combo_matrix="data/{dataset}/helper_model_{n_individuals}{subpopulation}_combo_matrix.npy",
        #tricky_variants="data/{dataset}/tricky_variants.npy",
    output:
        genotypes="data/{dataset}/kageMappedReadsN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment,[A-Za-z0-9_]+}.vcf.gz",
        probs="data/{dataset}/kageMappedReadsN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment,[A-Za-z0-9_]+}.vcf.gz.tmp.probs.npy",
        count_probs="data/{dataset}/kageMappedReadsN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment,[A-Za-z0-9_]+}.vcf.gz.tmp.count_probs.npy",
        benchmark_report="data/{dataset}/benchmarks/kageMappedReadsN{n_individuals}{subpopulation,[a-z]+}_{experiment}.tsv"
    threads:
        config["n_threads"]
    benchmark:
        "data/{dataset}/benchmarks/kageMappedReadsN{n_individuals,\d+}{subpopulation,[a-z]+}_{experiment}-snakemake.tsv"
    
    params:
        sample_name=get_sample_name_from_experiment,
        read_coverage=get_read_coverage_from_experiment
    shell:
        "/usr/bin/time -v kage genotype -c {input.node_counts} "
        "-g {input.variant_to_nodes} " 
        "-v {input.variants} " 
        "-A {input.model} " 
        "-f {input.helper_model} "
        "-F {input.helper_model_combo_matrix} "
        #"-x {input.tricky_variants} "
        "-C CombinationModelGenotyper "
        "-t {config[n_threads]} "
        "-q 0.5 "
        "--average-coverage {params.read_coverage} "
        "-o {output.genotypes}.tmp "
        "--sample-name-output {params.sample_name} 2> {output.benchmark_report} "
        "&& bgzip -c {output.genotypes}.tmp > {output.genotypes} "

