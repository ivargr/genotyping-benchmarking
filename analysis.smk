include:
    "kage.smk"

include:
    "genotype_methods.smk"

include:
    "pangenie.smk"

def get_truth_file_vcf_url(wildcards):
    return config["truth_datasets"][wildcards.truth_dataset]["vcf_url"]

def get_truth_file_regions_url(wildcards):
    return config["truth_datasets"][wildcards.truth_dataset]["regions_url"]

def get_dataset_regions_comma_separated(wildcards):

    return config["analysis_regions"]["dataset" + wildcards.n]["region"].replace(" ", ",")

def get_genome_size_simulated_data_set(wildcards):
    return config["analysis_regions"]["simulated_dataset" + wildcards.number]["genome_size"]


# hack to make bcftools --samples case inssensitive by providing the sample in upper and lower case
def samples_list(wildcards):
    return wildcards.truth_dataset.lower() + "," + wildcards.truth_dataset.upper()

ruleorder:
    debug_genotyping_mapped_reads > debug_genotyping


rule download_truth_file:
    output:
        vcf="data/{dataset}/original_{truth_dataset}.vcf.gz",
        vcf_index="data/{dataset}/original_{truth_dataset}.vcf.gz.tbi",
        regions_file="data/{dataset}/original_{truth_dataset}_regions.bed"

    params:
        vcf_url=get_truth_file_vcf_url,
        regions_url=get_truth_file_regions_url
    conda: "envs/prepare_data.yml"
    shell:
        # hacky way to convert chromosomes to numeric:
        "wget -O {output.vcf}.tmp.gz {params.vcf_url} && gunzip -f {output.vcf}.tmp.gz && sed -i 's/chr//g' {output.vcf}.tmp && bgzip -f -c {output.vcf}.tmp > {output.vcf} && tabix -f -p vcf {output.vcf} && "
        "wget -O {output.regions_file} {params.regions_url} && sed -i 's/chr//g' {output.regions_file}"



# whole region is in regions file
rule create_simulated_regions_file:
    output:
        regions_file="data/simulated_dataset{number}/truth_{truth_dataset}_regions.bed"
    params:
        genome_size=get_genome_size_simulated_data_set
    shell:
        "echo '1\t1\t{params.genome_size}' > {output}"


"""
rule create_simulated_truth_file:
    input:
        population_vcf="data/simulated_dataset{number}/variants.vcf"
    output:
        vcf="data/simulated_dataset{number}/truth_seed{seed,\d+}.vcf",
        vcf_gz="data/simulated_dataset{number}/truth_seed{seed,\d+}.vcf.gz",
    shell:
        "graph_read_simulator simulate_individual_vcf -s {wildcards.seed} -v {input} -o {output.vcf} && bgzip -f -c {output.vcf} > {output.vcf_gz} && tabix -f -p vcf {output.vcf_gz}"
"""




# makes a truth-file for non-simulated variants, i.e. variants downloaded
rule create_truth_file:
    input:
        vcf="data/dataset{n}/original_{truth_dataset}.vcf.gz",
        regions_file="data/dataset{n}/original_{truth_dataset}_regions.bed"

    output:
        vcf="data/dataset{n,\d+}/truth_{truth_dataset}.vcf.gz",
        regions_file="data/dataset{n,\d+}/truth_{truth_dataset}_regions.bed"
    params:
        regions = get_dataset_regions_comma_separated,
        samples=samples_list
    conda: "envs/prepare_data.yml"
    shell:
        "bcftools view --regions {params.regions} -f PASS --force-samples --samples {params.samples} {input.vcf} |  "
        "python3 scripts/filter_uncertain_sv.py | bgzip -c -f > {output.vcf} "
        "&& tabix -f -p vcf {output.vcf} && "
        "python scripts/extract_regions_from_bed.py {input.regions_file} {params.regions} > {output.regions_file}"


rule create_truth_file_without_long_indels:
    input: "data/{dataset}/truth_{truth_dataset}.vcf.gz"
    output: "data/{dataset}/truthShortIndels_{truth_dataset}.vcf.gz"
    conda: "envs/prepare_data.yml"
    shell: "zcat {input} | python scripts/filter_vcf_on_indel_length.py 20 | bgzip -c > {output}"


rule intersect_truth_file_with_callable_variants:
    input:
        truth="data/{dataset}/truth_{truth_dataset}.vcf.gz",
        vcf="data/{dataset}/variants_no_genotypes.vcf"
    output:
        "data/{dataset}/truthOnlyCallable_{truth_dataset}.vcf.gz"
    
    shell:
        "obgraph intersect_vcfs -a {input.truth} -b {input.vcf} -o {output}.tmp && bgzip -f -c {output}.tmp > {output}"


rule run_happy:
    input:
        genotypes="data/{dataset}/{run}.vcf.gz",
        truth_vcf="data/{dataset}/truth_{truth_dataset}.vcf.gz",
        #truth_vcf_only_callable="data/{dataset}/truthOnlyCallable_{truth_dataset}.vcf.gz",
        #truth_short_indels="data/{dataset}/truthShortIndels_{truth_dataset}.vcf.gz",
        truth_regions_file="data/{dataset}/truth_{truth_dataset}_regions.bed",
        ref="data/{dataset}/ref.fa"
    output:
        output_file="data/{dataset}/happy-{truth_dataset}-{run}.extended.csv",
        summary_output_file="data/{dataset}/happy-{truth_dataset}-{run}.summary.csv",
        #output_file_only_callable="data/{dataset}/happy-{truth_dataset}-{run}-only-callable.extended.csv",
        #output_file_short_indels="data/{dataset}/happy-{truth_dataset}-{run}-short-indels.extended.csv",
    conda:
        "envs/happy.yml"
    shell:
        """zcat {input.genotypes} | awk '{{ $6 = ($6 == "inf" ? 0 : $6) }} 1' OFS="\\t" | bgzip -c > {input.genotypes}.inf-replaced.vcf.gz && """
        "tabix -f -p vcf {input.genotypes}.inf-replaced.vcf.gz && "
        # skip long indels
        #"hap.py {input.truth_short_indels} {input.genotypes}.inf-replaced.vcf.gz --no-leftshift -r {input.ref} -o data/{wildcards.dataset}/happy-{wildcards.truth_dataset}-{wildcards.run}-short-indels -f {input.truth_regions_file} --no-decompose --engine=vcfeval && "
        # normal
        "hap.py {input.truth_vcf} {input.genotypes}.inf-replaced.vcf.gz --no-leftshift -r {input.ref} -o data/{wildcards.dataset}/happy-{wildcards.truth_dataset}-{wildcards.run} -f {input.truth_regions_file} --no-decompose --engine=vcfeval "
        # check against only truth variants that are in input variant set:
        #"hap.py {input.truth_vcf_only_callable} {input.genotypes}.inf-replaced.vcf.gz --no-leftshift -r {input.ref} -o data/{wildcards.dataset}/happy-{wildcards.truth_dataset}-{wildcards.run}-only-callable -f {input.truth_regions_file} --no-decompose --engine=vcfeval "
        #"hap.py {input.truth_vcf} {input.genotypes}.inf-replaced.vcf.gz -r {input.ref} -o data/{wildcards.dataset}/happy-{wildcards.truth_dataset}-{wildcards.run} -f {input.truth_regions_file}"

"""
rule bgzip_result_file:
    input: "data/{dataset}/{method}_{reads}.{n_individuals,\d+}individuals.vcf"
    output: "data/{dataset}/{method}_{reads}.{n_individuals,\d+}individuals.vcf.gz"
    shell: "bgzip -f -c {input} > {output} && tabix -f -p vcf {output}"
"""


rule add_sampe_name_to_vcf_for_trio_analysis:
    input:
        "data/{dataset}/{method}_{sample_name}_real_reads_{coverage,\d+}x.vcf.gz"
    output:
        "data/{dataset}/{method,[a-zA-Z0-9]+}_{sample_name}_real_reads_{coverage,\d+}x.with_sample_name.vcf.gz"
    conda: "envs/prepare_data.yml"
    shell:
        "echo '{wildcards.sample_name}' > {wildcards.method}_{wildcards.sample_name}.sample_name && "
        "bcftools reheader --samples {wildcards.method}_{wildcards.sample_name}.sample_name -o {output} {input} && "
        "tabix -f -p vcf {output}"

rule analyse_trio_concordance:
    input:
        child="data/{dataset}/{method}_hg002_real_reads_15x.with_sample_name.vcf.gz",
        mother="data/{dataset}/{method}_hg004_real_reads_15x.with_sample_name.vcf.gz",
        father="data/{dataset}/{method}_hg003_real_reads_15x.with_sample_name.vcf.gz",
        reference="data/{dataset}/ref.fa"
    output:
        pedigree_file="data/{dataset}/aj_pedigree_{method}/aj_pedigree_tab_delim_detailed_log.tsv"
    params:
        outputdir="data/{dataset}/aj_pedigree_{method}"
    shell:
        """
        mkdir -p {params.outputdir} && 
        export LD_LIBRARY_PATH=VBT-TrioAnalysis/lib/ && 
        vbt mendelian -ref {input.reference} -mother {input.mother} -father {input.father} \
        -child {input.child} -pedigree resources/pedigree_aj.ped -outDir {params.outputdir} -out-prefix aj_pedigree --output-violation-regions
        """


rule debug_genotyping:
    input:
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        kmer_index="data/{dataset}/kmer_index_only_variants.npz",
        reverse_variant_kmers="data/{dataset}/reverse_variant_kmers.npz",
        input_vcf="data/{dataset}/variants_no_genotypes.vcf",
        genotyped_vcf="data/{dataset}/{method}N{n_individuals,\d+}all_{experiment}.vcf.gz",
        truth_vcf="data/{dataset}/truth_{truth_id}.vcf.gz",
        truth_regions="data/{dataset}/truth_{truth_id}_regions.bed",
        node_counts="data/{dataset}/{experiment}.I1000.node_counts.npy",
        #model="data/{dataset}/combination_model.npz",
        model="data/{dataset}/sampling_count_model_{n_individuals}all.npz",
        helper_variants="data/{dataset}/helper_model_{n_individuals}all.npy",
        combination_matrix="data/{dataset}/helper_model_{n_individuals}all_combo_matrix.npy",
        genotype_probs="data/{dataset}/{method}N{n_individuals}all_{experiment}.vcf.gz.tmp.probs.npy",
        genotype_count_probs="data/{dataset}/{method}N{n_individuals}all_{experiment}.vcf.gz.tmp.count_probs.npy",
        #pangenie="data/{dataset}/pangenieN32_{experiment}.vcf.gz"

    output:
        "data/{dataset}/debugging-{method}N{n_individuals,\d+}all-{truth_id}-{experiment}.txt"

    shell:
        "kage analyse_variants -g {input.variant_to_nodes} "
        "-i {input.kmer_index} "
        "-R {input.reverse_variant_kmers} "
        "-k 31 "
        "-v {input.input_vcf} "
        "-P {input.genotyped_vcf} "
        "-T {input.truth_vcf} "
        "-n {input.node_counts} "
        "-m {input.model} "
        "-f {input.helper_variants} "
        "-F {input.combination_matrix} "
        "-p {input.genotype_probs} "
        "-c {input.genotype_count_probs} "
        #"-a {input.pangenie} "
        "-t {input.truth_regions} > {output}"



rule debug_genotyping_mapped_reads:
    input:
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        kmer_index="data/{dataset}/kmer_index_only_variants.npz",
        reverse_variant_kmers="data/{dataset}/reverse_variant_kmers.npz",
        input_vcf="data/{dataset}/variants_no_genotypes.vcf",
        genotyped_vcf="data/{dataset}/kageMappedReadsN{n_individuals,\d+}all_{experiment}.vcf.gz",
        truth_vcf="data/{dataset}/truth_{truth_id}.vcf.gz",
        truth_regions="data/{dataset}/truth_{truth_id}_regions.bed",
        node_counts="data/{dataset}/{experiment}.node_counts_giraffe.npy",
        model="data/{dataset}/model_for_read_mapping_with_sampled_reads.npz",
        helper_variants="data/{dataset}/helper_model_{n_individuals}all.npy",
        combination_matrix="data/{dataset}/helper_model_{n_individuals}all_combo_matrix.npy",
        genotype_probs="data/{dataset}/kageMappedReadsN{n_individuals}all_{experiment}.vcf.gz.tmp.probs.npy",
        genotype_count_probs="data/{dataset}/kageMappedReadsN{n_individuals}all_{experiment}.vcf.gz.tmp.count_probs.npy",
        #pangenie="data/{dataset}/pangenieN32_{experiment}.vcf.gz"

    output:
        "data/{dataset}/debugging-kageMappedReadsN{n_individuals,\d+}all-{truth_id}-{experiment}.txt"

    shell:
        "kage analyse_variants -g {input.variant_to_nodes} "
        "-i {input.kmer_index} "
        "-R {input.reverse_variant_kmers} "
        "-k 31 "
        "-v {input.input_vcf} "
        "-P {input.genotyped_vcf} "
        "-T {input.truth_vcf} "
        "-n {input.node_counts} "
        "-m {input.model} "
        "-f {input.helper_variants} "
        "-F {input.combination_matrix} "
        "-p {input.genotype_probs} "
        "-c {input.genotype_count_probs} "
        #"-a {input.pangenie} "
        "-t {input.truth_regions} > {output}"


rule find_variants_with_nonunique_kmers2:
    input:
        variants="data/{dataset}/variants_no_genotypes.vcf",
        variant_to_nodes="data/{dataset}/variant_to_nodes.npz",
        reverse_kmer_index="data/{dataset}/reverse_variant_kmers.npz",
        population_kmers="data/{dataset}/kmer_index.npz"
    output:
        "data/{dataset}/variants_with_nonunique_kmers.pickle"
    
    shell:
        "python3 scripts/find_variants_with_nonunique_kmers.py "
        "{input.variants} "
        "{input.reverse_kmer_index} "
        "{input.population_kmers} "
        "{input.variant_to_nodes} "
        "{output} "


rule analyse_variants_with_nonunique_kmers:
    input:
        variants="data/{dataset}/variants_no_genotypes.vcf",
        naive_kage="data/{dataset}/naivekageN2548all_hg002_simulated_reads_15x.vcf.gz",
        kage_no_helper_model="data/{dataset}/kageNoHelperModelN250all_hg002_simulated_reads_15x.vcf.gz",
        kage="data/{dataset}/usN2548all_hg002_simulated_reads_15x.vcf.gz",
        truth="data/{dataset}/truth_hg002.vcf.gz",
        nonunique="data/{dataset}/variants_with_nonunique_kmers.pickle",
       
    output:
        "data/{dataset}/nonunique_variants_report.txt"
        
    
    shell:
        "python3 scripts/analyse_variants_with_nonunnique_kmers.py "
        "{input.variants} "
        "{input.naive_kage} "
        "{input.kage_no_helper_model} "
        "{input.kage} "
        "{input.truth} "
        "{input.nonunique} > {output}"
