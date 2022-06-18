
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
        "obgraph merge_graphs -o {output} -g {input}"


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


rule make_haplotype_matrix:
    input:
        vcf="data/{dataset}/variants_{n_individuals}{subpopulation}.vcf.gz",
        vcf_no_genotypes="data/{dataset}/variants_no_genotypes.vcf"
    output:
        haplotype_matrix="data/{dataset}/haplotype_matrix_{n_individuals,\d+}{subpopulation,[a-z]+}.npy",
    threads: config["n_threads_data_quarter"]
    benchmark: "data/{dataset}/benchmarks/make_haplotype_matrix_{n_individuals}{subpopulation}.tsv"
    resources:
        mem_gb=300
    conda: "envs/kage.yml"
    shell:
        "n_variants=$(grep -v '#' {input.vcf_no_genotypes} | wc -l) || true && "
        "obgraph make_haplotype_matrix  -v {input.vcf} -o {output.haplotype_matrix} -n {wildcards.n_individuals} -m $n_variants -t {config[n_threads_data]}"


rule make_genotype_matrix:
    input:
        vcf="data/{dataset}/variants_{n_individuals}{subpopulation}.vcf.gz",
        vcf_no_genotypes="data/{dataset}/variants_no_genotypes.vcf"
    output:
        genotype_matrix="data/{dataset}/genotype_matrix_{n_individuals,\d+}{subpopulation,[a-z]+}.npy",
    threads: config["n_threads_data_quarter"]
    benchmark: "data/{dataset}/benchmarks/make_genotype_matrix_{n_individuals}{subpopulation}.tsv"
    resources:
        mem_gb=300
    conda: "envs/kage.yml"
    shell:
        #"n_individuals=$(zcat {input.vcf} | head -n 1000 | grep -i '#chrom' | python3 -c 'import sys; print(len(list(sys.stdin)[0].split())-9)') || true && "
        "n_variants=$(grep -v '#' {input.vcf_no_genotypes} | wc -l) || true && "
        "obgraph make_genotype_matrix -v {input.vcf} -o {output.genotype_matrix} -n {wildcards.n_individuals} -m $n_variants -t {config[n_threads_data]}"
    
    
rule make_genotype_txt_matrix:
    input:
        vcf="data/{dataset}/variants_{n_individuals}{subpopulation}.vcf.gz",
    output:
        genotype_matrix="data/{dataset}/genotype_matrix_{n_individuals,\d+}{subpopulation,[a-z]+}.txt.gz",
    threads: 1
    benchmark: "data/{dataset}/benchmarks/make_genotype_txt_matrix_{n_individuals}{subpopulation}.tsv"
    conda:
        "envs/prepare_data.yml"
    shell:
        'gunzip -c {input} | grep -v "^#" | cut -f 10- | tr -d "|" | tr -d "\n" | tr -d "\t" | pigz -c > {output} '


rule make_phased_genotype_matrix:
    input:
        #vcf="data/{dataset}/variants_{n_individuals}{subpopulation}.vcf.gz",
        vcf="data/{dataset}/genotype_matrix_{n_individuals}{subpopulation}.txt.gz",
        vcf_no_genotypes="data/{dataset}/variants_no_genotypes.vcf"
    output:
        genotype_matrix="data/{dataset}/phased_genotype_matrix_{n_individuals,\d+}{subpopulation,[a-z]+}.npz",
    threads: 1
    benchmark: "data/{dataset}/benchmarks/make_phased_genotype_matrix_{n_individuals}{subpopulation}.tsv"
    resources:
        mem_gb=150
    conda: "envs/kage.yml"
    shell:
        #"n_individuals=$(zcat {input.vcf} | head -n 1000 | grep -i '#chrom' | python3 -c 'import sys; print(len(list(sys.stdin)[0].split())-9)') || true && "
        "n_variants=$(grep -v '#' {input.vcf_no_genotypes} | wc -l) || true && "
        "obgraph make_genotype_matrix --make-phased-matrix True -v {input.vcf} -o {output.genotype_matrix} -n {wildcards.n_individuals} -m $n_variants -t {config[n_threads_data]}"


rule make_kmer_counter_index:
    input:
        kmer_index="data/{dataset}/kmer_index_only_variants_with_revcomp.npz"
    output:
        index="data/{dataset}/counter_index_only_variants_with_revcomp.npz"
    shell:
        "graph_kmer_index create_counter_index -i {input.kmer_index} -o {output.index}"


rule make_count_model:
    input:
        haplotype_to_nodes="data/{dataset}/disc_backed_haplotype_to_nodes_{n_individuals}{subpopulation}.npz",
        graph="data/{dataset}/obgraph.npz",
        #counter_index="data/{dataset}/counter_index_only_variants_with_revcomp.npz"
        counter_index="data/{dataset}/kmer_index_only_variants_with_revcomp.npz"
    output:
        "data/{dataset}/sampling_count_model_{n_individuals,\d+}{subpopulation}.npz"
    threads:
        config["n_threads_data_quarter"]
    resources:
        mem_gb=450
    benchmark: "data/{dataset}/benchmarks/make_count_model_{n_individuals}{subpopulation}.tsv"
    shell:
        "kage sample_node_counts_from_population -g {input.graph} "
        "-H {input.haplotype_to_nodes} "
        #"-H data/dataset1/disc_backed_haplotype_to_nodes "
        "-i {input.counter_index} -o {output} -t 25 --max-count 15"


rule make_haplotype_to_nodes:
    input:
        #graph="data/{dataset}/obgraph.npz",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        #vcf = "data/{dataset}/variants_{n_individuals}{subpopulation}.vcf.gz",
        phased_genotype_matrix="data/{dataset}/phased_genotype_matrix_{n_individuals}{subpopulation}.npz"
    output:
        "data/{dataset}/haplotype_to_nodes_{n_individuals,\d+}{subpopulation}.npz"
    benchmark:
        "data/{dataset}/benchmarks/make_haplotype_to_nodes_{n_individuals}{subpopulation}.tsv"
    shell:
        "n_haplotypes=$((2*{wildcards.n_individuals})) && "
        #"obgraph make_haplotype_to_nodes -g {input.graph} -v {input.vcf} -n $n_haplotypes -o {output}"
        "obgraph make_haplotype_to_nodes_bnp -g {input.variant_to_nodes} -v {input.phased_genotype_matrix} -n $n_haplotypes -o {output}"


rule make_disc_backed_haplotype_to_nodes:
    input:
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        phased_genotype_matrix="data/{dataset}/phased_genotype_matrix_{n_individuals}{subpopulation}.npz"
    output:
        "data/{dataset}/disc_backed_haplotype_to_nodes_{n_individuals,\d+}{subpopulation}.npz"
    benchmark:
        "data/{dataset}/benchmarks/make_haplotype_to_nodes_{n_individuals}{subpopulation}.tsv"
    shell:
        "n_haplotypes=$((2*{wildcards.n_individuals})) && "
        #"obgraph make_haplotype_to_nodes -g {input.graph} -v {input.vcf} -n $n_haplotypes -o {output}"
        "obgraph make_haplotype_to_nodes_bnp -g {input.variant_to_nodes} -v {input.phased_genotype_matrix} -n $n_haplotypes -o {output} -d True"


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
        genotype_matrix="data/{dataset}/genotype_matrix_{n_individuals}{subpopulation}.npy",
        most_similar_variant_lookup="data/{dataset}/most_similar_variant_lookup_{n_individuals}{subpopulation}.npz",
    output:
        transition_probs = "data/{dataset}/transition_probs_{n_individuals,\d+}{subpopulation,[a-z]+}.npy",
    threads: config["n_threads_data_quarter"]
    benchmark: "data/{dataset}/benchmarks/make_transition_probabilities_{n_individuals}{subpopulation}.tsv"
    resources: mem_gb=100
    conda: "envs/kage.yml"
    shell:
        "obgraph make_transition_probabilities -t {config[n_threads_data_quarter]} -G {input.genotype_matrix} -m {input.most_similar_variant_lookup} -o {output.transition_probs}"


rule make_most_similar_variant_lookup:
    input:
        genotype_matrix="data/{dataset}/genotype_matrix_{n_individuals}{subpopulation}.npy",
    output:
        most_similar_variant_lookup="data/{dataset}/most_similar_variant_lookup_{n_individuals,\d+}{subpopulation,[a-z]+}.npz",
    threads: config["n_threads_data_quarter"]
    benchmark: "data/{dataset}/benchmarks/make_most_similar_variant_lookup_{n_individuals}{subpopulation}.tsv"
    resources: mem_gb=100
    conda: "envs/kage.yml"
    shell:
        "obgraph analyse_genotype_matrix -t {config[n_threads_data_quarter]} -G {input.genotype_matrix} -o {output.most_similar_variant_lookup}"


rule get_genotype_frequencies:
    input:
        genotype_matrix = "data/{dataset}/genotype_matrix_{n_individuals}{subpopulation}.npy",
    output:
        genotype_frequencies="data/{dataset}/genotype_frequencies_{n_individuals,\d+}{subpopulation,[a-z]+}.npz"
    threads: config["n_threads_data"]
    benchmark: "data/{dataset}/benchmarks/get_genotype_frequencies_{n_individuals}{subpopulation}.tsv"
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
    params:
        genome_size=get_dataset_genome_size
    resources:
        mem_gb=300
    threads: config["n_threads_data_quarter"]
    benchmark:
        "data/{dataset}/benchmarks/sample_kmers_from_linear_reference.tsv"
    conda: "envs/kage.yml"
    shell:
        "graph_kmer_index make -s 1 -k {config[k]} --include-reverse-complement True -G {params.genome_size} -o data/{wildcards.dataset}/linear_kmers -R {input} -n ref -t {config[n_threads_data_quarter]} "


rule make_linear_reference_kmer_index:
    input:
        flat = "data/{dataset}/linear_kmers.npz",
    output:
        index="data/{dataset}/linear_kmer_index.npz"
    conda:
        "envs/kage.yml"
    shell:
        "graph_kmer_index make_from_flat -o {output.index} -f {input.flat} -m 20000033"


rule make_linear_reference_kmer_index_counter:
    input:
        flat="data/{dataset}/linear_kmers.npz",
    output:
        index="data/{dataset}/linear_kmers_counter.npz"
    benchmark: "data/{dataset}/make_linear_reference_kmer_index_counter.tsv"
    conda:
        "envs/kage.yml"
    shell:
        "graph_kmer_index count_kmers -f {input.flat} -o {output.index} --subsample-ratio 1"


rule add_reverse_complements_to_kmers:
    input: "data/{dataset}/{kmers}.npz"
    output: "data/{dataset}/{kmers}_with_reverse_complements.npz"
    conda: "envs/kage.yml"
    shell: "graph_kmer_index add_reverse_complements -f {input} -k {config[k]} -o {output}"


rule make_helper_model:
    input:
        genotype_matrix="data/{dataset}/genotype_matrix_{n_individuals}{subpopulation}.npy",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        node_count_model="data/{dataset}/combination_model.npz",
    output:
        helper_model="data/{dataset}/helper_model_{n_individuals,\d+}{subpopulation,[a-z]+}.npy",
        helper_model_combo_matrix="data/{dataset}/helper_model_{n_individuals,\d+}{subpopulation,[a-z]+}_combo_matrix.npy"
    threads: config["n_threads_data"]
    resources: mem_gb=300
    conda: "envs/kage.yml"
    shell:
        "kage create_helper_model -o data/{wildcards.dataset}/helper_model_{wildcards.n_individuals}{wildcards.subpopulation} "
        "-g {input.genotype_matrix} -w 100 -v {input.variant_to_nodes} -n {input.node_count_model} -t {config[n_threads_data_quarter]} "

rule make_critical_paths_index:
    input:
        graph = "data/{dataset}/obgraph.npz",
    output:
        index = "data/{dataset}/critical_paths2.npz",
    benchmark: "data/{dataset}/make_critical_paths_index.tsv"
    shell:
        "graph_kmer_index find_critical_paths -g {input.graph} -k {config[k]} -o {output.index}"



rule make_position_id_index:
    input:
        graph = "data/{dataset}/obgraph.npz",
    output:
        index = "data/{dataset}/position_id_index.npz",
    benchmark: "data/{dataset}/make_position_id_index.tsv"
    shell:
        "obgraph make_position_id -g {input.graph} -o {output.index}"


rule get_variant_kmers:
    input:
        graph="data/{dataset}/obgraph.npz",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        vcf="data/{dataset}/variants_no_genotypes.vcf",
        #linear_kmer_index="data/{dataset}/linear_kmer_index.npz",
        linear_kmer_index="data/{dataset}/linear_kmers_counter.npz",
        position_id_index="data/{dataset}/position_id_index.npz"
    output:
        variant_kmers="data/{dataset}/variant_kmers_short_variants.npz",
    benchmark:
        "data/{dataset}/benchmarks/get_variant_kmers.tsv"
    resources: mem_gb=400
    threads: config["n_threads_data"]
    conda: "envs/kage.yml"
    shell:
        "graph_kmer_index make_unique_variant_kmers -g {input.graph} -V {input.variant_to_nodes} -k {config[k]} -o {output} -v {input.vcf} "
        " -t {config[n_threads_data_quarter]} -c 20000 --max-variant-nodes 3 -I {input.linear_kmer_index} -p {input.position_id_index} -D True "

rule get_structural_variant_kmers:
    input:
        graph="data/{dataset}/obgraph.npz",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        vcf="data/{dataset}/variants_no_genotypes.vcf",
        #linear_kmer_index="data/{dataset}/linear_kmer_index.npz",
        linear_kmer_index="data/{dataset}/linear_kmers_counter.npz",
        position_id_index="data/{dataset}/position_id_index.npz"
    output:
        variant_kmers="data/{dataset}/structural_variant_kmers.npz",
    benchmark:
        "data/{dataset}/benchmarks/make_variant_kmer_index.tsv"
    resources: mem_gb=400
    threads: config["n_threads_data"]
    conda: "envs/kage.yml"
    shell:
        "graph_kmer_index sample_kmers_from_structural_variants -g {input.graph} -V {input.variant_to_nodes} -i {input.linear_kmer_index} -k {config[k]} -o {output}"


rule merge_structural_variant_kmers_with_variant_kmers:
    input:
        variant_kmers="data/{dataset}/variant_kmers_short_variants.npz",
        structural_variant_kmers="data/{dataset}/structural_variant_kmers.npz",
    output:
        index="data/{dataset}/variant_kmers.npz",
    benchmark:
        "data/{dataset}/benchmarks/merge_structural_variant_kmers_with_variant_kmers.tsv"
    resources: mem_gb=50
    threads: 1
    conda: "envs/kage.yml"
    shell:
        "graph_kmer_index merge_flat_kmers --flat-kmers {input.variant_kmers},{input.structural_variant_kmers} -o {output}"


rule make_reverse_variant_kmer_index:
    input:
        variant_kmers="data/{dataset}/variant_kmers.npz",
    output:
        reverse_index="data/{dataset}/reverse_variant_kmers.npz"
    benchmark: "data/{dataset}/benchmarks/make_reverse_variant_kmer_index.tsv"
    shell:
        "graph_kmer_index make_reverse -f {input} -o {output}"


rule make_variant_kmer_index_with_reverse_complements:
    input:
        variant_kmers = "data/{dataset}/variant_kmers.npz",
    output:
        index = "data/{dataset}/kmer_index_only_variants_with_revcomp.npz",
    shell:
        "graph_kmer_index make_from_flat -o {output.index} -f {input.variant_kmers} -m 200000033 -k {config[k]} -r True"


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
        #variant_kmers="data/{dataset}/variant_kmers.npz",
        variant_kmers="data/{dataset}/kmer_index_only_variants_with_revcomp.npz",
        position_id_index="data/{dataset}/position_id_index.npz",
        critical_paths="data/{dataset}/critical_paths2.npz"
    params:
        genome_size=get_dataset_genome_size
    output:
        flat_kmers="data/{dataset}/flat_kmers_for_model.npz",
    benchmark: "data/{dataset}/benchmarks/sample_kmers_for_kmer_model.tsv"
    resources: mem_gb=100
    threads: config["n_threads_data"]
    conda: "envs/kage.yml"
    shell:
        #"graph_kmer_index make -s 1 -k {config[k]} -G {params.genome_size} -o {output.flat_kmers} -g {input.graph} --max-frequency 100000 -v 5 -t {config[n_threads_data]} --only-save-one-node-per-kmer True --include-reverse-complement True --whitelist {input.variant_kmers}"
        "graph_kmer_index index -k {config[k]} -g {input.graph} -t 25 --whitelist {input.variant_kmers} "
        " -r True -p {input.position_id_index} -c {input.critical_paths} -o {output.flat_kmers} --max-variant-nodes 3"


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
        "kage model -i {input.kmer_index} -g {input.variant_to_nodes} -r {input.reverse_variant_kmers} -m {params.n_nodes} -o {output.haplotype_model} -k {config[k]} -t {config[n_threads_data_quarter]} -V v3 && "
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
        #model= "data/{dataset}/combination_model.npz",
        model= "data/{dataset}/sampling_count_model_{n_individuals}{subpopulation}.npz",
        reverse_variant_kmers="data/{dataset}/reverse_variant_kmers.npz",
    output:
        tricky_variants="data/{dataset}/tricky_variants_{n_individuals,\d+}{subpopulation}.npy",
        tricky_variants_nonunique_kmers="data/{dataset}/tricky_variants_{n_individuals}{subpopulation}_nonunique_kmers.npy"
    conda: "envs/kage.yml"
    shell:
        "kage find_tricky_variants -v {input.variant_to_nodes} -m {input.model} -r {input.reverse_variant_kmers} -o {output.tricky_variants} -M 1000 && "
        "kage find_tricky_variants -v {input.variant_to_nodes} -m {input.model} -r {input.reverse_variant_kmers} -o {output.tricky_variants_nonunique_kmers} -u True"


rule find_variants_with_nonunique_kmers:
    input:
        kmer_index="data/{dataset}/kmer_index.npz",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        reverse_variant_kmers="data/{dataset}/reverse_variant_kmers.npz",
    output:
        tricky_variants="data/{dataset}/variants_with_nonunique_kmers.npy",
    conda: "envs/kage.yml"
    shell:
        "kage find_variants_with_nonunique_kmers -v {input.variant_to_nodes} -i {input.kmer_index} -r {input.reverse_variant_kmers} -o {output.tricky_variants} "


rule make_index_bundle:
    input:
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        numpy_variants="data/{dataset}/numpy_variants.npz",
        model="data/{dataset}/combination_model.npz",
        tricky_variants="data/{dataset}/tricky_variants.npy",
        helper_variants="data/{dataset}/helper_model_{n_individuals}{subpopulation}.npy",
        combo_matrix="data/{dataset}/helper_model_{n_individuals}{subpopulation}_combo_matrix.npy",
        kmer_index="data/{dataset}/kmer_index_only_variants_with_revcomp.npz"
    output:
        "data/{dataset}/index_{n_individuals,\d+}{subpopulation,[a-z]+}.npz"
    conda: "envs/kage.yml"
    shell: 
        "kage make_index_bundle "
        "-g {input.variant_to_nodes} "
        "-v {input.numpy_variants} "
        "-A {input.model} "
        "-x {input.tricky_variants} "
        "-f {input.helper_variants} "
        "-F {input.combo_matrix} "
        "-i {input.kmer_index} "
        "-o {output} "


