include:
    "alignment_free_graph_genotyper.smk"

include:
    "genotype_methods.smk"

def get_truth_file_vcf_url(wildcards):
    return config["truth_datasets"][wildcards.truth_dataset]["vcf_url"]

def get_truth_file_regions_url(wildcards):
    return config["truth_datasets"][wildcards.truth_dataset]["regions_url"]

def get_dataset_regions_comma_separated(wildcards):
    return config["analysis_regions"][wildcards.dataset]["region"].replace(" ", ",")

rule download_truth_file:
    output:
        vcf="data/{dataset}/original_{truth_dataset}.vcf.gz",
        vcf_index="data/{dataset}/original_{truth_dataset}.vcf.gz.tbi",
        regions_file="data/{dataset}/original_{truth_dataset}_regions.bed"

    params:
        vcf_url=get_truth_file_vcf_url,
        regions_url=get_truth_file_regions_url

    shell:
        # hacky way to convert chromosomes to numeric:
        "wget -O {output.vcf}.tmp.gz {params.vcf_url} && gunzip -f {output.vcf}.tmp.gz && sed -i 's/chr//g' {output.vcf}.tmp && bgzip -f -c {output.vcf}.tmp > {output.vcf} && tabix -p vcf {output.vcf} && "
        "wget -O {output.regions_file} {params.regions_url} && sed -i 's/chr//g' {output.regions_file}"


rule create_truth_file:
    input:
        vcf="data/{dataset}/original_{truth_dataset}.vcf.gz",
        regions_file="data/{dataset}/original_{truth_dataset}_regions.bed"

    output:
        vcf="data/{dataset}/truth_{truth_dataset}.vcf.gz",
        regions_file="data/{dataset}/truth_{truth_dataset}_regions.bed"

    params:
        regions = get_dataset_regions_comma_separated

    shell:
        "bcftools view -O z --regions {params.regions} {input.vcf} > {output.vcf} && tabix -p vcf {output.vcf} && "
        "python3 scripts/extract_regions_from_bed.py {input.regions_file} {params.regions} > {output.regions_file}"


rule run_happy:
    input:
        genotypes="data/{dataset}/{run}.vcf.gz",
        truth_vcf="data/{dataset}/truth_{truth_dataset}.vcf.gz",
        truth_regions_file="data/{dataset}/truth_{truth_dataset}_regions.bed",
        ref="data/{dataset}/ref.fa"

    output:
        output_file="data/{dataset}/happy-{truth_dataset}-{run}.summary.csv"

    shell:
        "/root/hap.py-install/bin/hap.py {input.truth_vcf} {input.genotypes} --no-leftshift -r {input.ref} -o data/{wildcards.dataset}/happy-{wildcards.truth_dataset}-{wildcards.run} -f {input.truth_regions_file} --no-decompose --engine=vcfeval"

"""
rule bgzip_result_file:
    input: "data/{dataset}/{method}_{reads}.{n_individuals,\d+}individuals.vcf"
    output: "data/{dataset}/{method}_{reads}.{n_individuals,\d+}individuals.vcf.gz"
    shell: "bgzip -f -c {input} > {output} && tabix -p vcf {output}"
"""






