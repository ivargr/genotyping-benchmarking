
def get_vg_graphs_name(wildcards):
    return " ".join(chr + ".vg" for chr in config["analysis_regions"]["dataset" + wildcards.number]["chromosomes"].split())


def get_vg_chromosome_graph_names(wildcards):
    return ["data/" + wildcards.dataset + "/graph_%s.vg"  % chromosome for chromosome in config["analysis_regions"][wildcards.dataset]["chromosomes"].split()]


rule make_vg_chromosome_graph:
    input:
        vcf="data/{dataset}/variants.vcf.gz",
        ref="data/{dataset}/ref.fa"
    output:
        graph="data/{dataset}/graph_{chromosome}.vg"
    shell:
        "vg construct --node-max 1000 -C -R {wildcards.chromosome} -a -r {input.ref} -v {input.vcf} > {output.graph}"


# also converts node id space
rule make_vg_xg_index:
    input:
        get_vg_chromosome_graph_names
    output:
        xg="data/{dataset}/graph.xg"
    shell:
        "vg ids -j {input} && vg index -L -t {config[n_threads]} -x {output} {input}"


rule make_chromosome_xg:
    input:
        graph="data/{dataset}/graph_{chromosome}.vg",
        whole_genome_xg="data/{dataset}/graph.xg"  # required this to make sure that node id conversion has been done
    output:
        xg="data/{dataset}/graph_{chromosome}.xg"
    shell:
        "vg index -x {output.xg} {input.graph}"


rule convert_chromosome_graph_to_gfa:
    input:
        xg="data/{dataset}/graph_{chromosome}.xg"
    output:
        gfa="data/{dataset}/graph_{chromosome}.gfa"
    shell:
        "vg view -g {input.xg} > {output.gfa}"


rule convert_whole_graph_to_gfa:
    input:
        "data/{dataset}/graph.xg"
    output:
        "data/{dataset}/graph.gfa"
    shell:
        "vg view -g {input} > {output}"

"""
rule make_obgraph_all_chromosomes:
    input:
        gfa="data/{dataset}/graph.gfa"
    output:
        "data/{dataset}/obgraph_no_dummy_nodes.npz"
    
    shell:
        "obgraph from_gfa -g {input.gfa} -o {output}"
"""

"""
rule add_indel_nodes:
    input:
        graph="data/{dataset}/obgraph_no_dummy_nodes.npz",
        variants="data/{dataset}/variants_no_genotypes.vcf"
    output:
        obgraph="data/{dataset}/obgraph.npz",
        edge_mapping="data/{dataset}/obgraph.npz.edge_mapping"
    
    shell:
        "obgraph add_indel_nodes -g {input.graph} -v {input.variants} -o {output.obgraph} && "
        "obgraph add_allele_frequencies -v {input.variants} -g {output.obgraph}"
"""
"""
rule make_obgraph_for_chromosome:
    input:
        xg="data/{dataset}/graph.xg",
        gfa="data/{dataset}/graph_{chromosome}.gfa"
    output:
        "data/{dataset}/obgraph_no_dummy_nodes_chr{chromosome}.npz"
    
    shell:
        "obgraph from_gfa -g {input.gfa} -o {output}"


rule add_dummy_nodes_to_obgraph:
    input:
        graph="data/{dataset}/obgraph_no_dummy_nodes_chr{chromosome}.npz",
        variants="data/{dataset}/variants_no_genotypes.vcf"
    output:
        "data/{dataset}/obgraph_chr{chromosome}.npz"
    
    shell:
        "obgraph add_indel_nodes -g {input.graph} -v {intput.variants} -o {output}"
        "obgraph add_allele_frequencies -c {wildcards.chromosome} -v {input.variants} -g {output}"
"""


rule make_vg_gbwt_index:
    input:
        xg="data/{dataset}/graph.xg",
        variants="data/{dataset}/variants.vcf.gz"
    output:
        "data/{dataset}/graph.gbwt"
    shell:
        "vg gbwt -x {input.xg} -o {output} --preset 1000gp --num-jobs {config[n_threads]} -v {input.variants}"


rule downsample_gbwt:
    input:
        xg="data/{dataset}/graph.xg",
        gbwt="data/{dataset}/graph.gbwt"
    output:
        gbwt="data/{dataset}/graph_downsampled.gbwt"
    shell:
        "vg gbwt -x {input.xg} -o {output.gbwt} -l {input.gbwt} -n 64"


rule make_gbwt_graph:
    input:
        xg="data/{dataset}/graph.xg",
        downsampled_gbwt="data/{dataset}/graph_downsampled.gbwt"
    output:
        gg="data/{dataset}/graph_downsampled.gg"
    shell:
        "vg gbwt -x {input.xg} -g {output.gg} {input.downsampled_gbwt}"


rule find_snarls:
    input:
        "data/{dataset}/graph.xg"
    output:
        "data/{dataset}/snarls.pb"
    shell:
        "vg snarls -T {input} > {output}"


rule make_dist_index:
    input:
        snarls="data/{dataset}/snarls.pb",
        xg="data/{dataset}/graph.xg"
    output:
        "data/{dataset}/graph.dist"
    shell:
        "vg index -t {config[n_threads]} -j {output} -s {input.snarls} {input.xg}"


rule make_minimizer_index:
    input:
        dist="data/{dataset}/graph.dist",
        gbwt="data/{dataset}/graph_downsampled.gbwt",
        gg="data/{dataset}/graph_downsampled.gg"
    output:
        "data/{dataset}/graph.min"
    shell:
        "vg minimizer -t {config[n_threads]} -p -i {output} -d {input.dist} -g {input.gbwt} -G {input.gg} "





