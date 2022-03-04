
include:
    "reference_genome.smk"

def get_dataset_regions(wildcards):
    return config["analysis_regions"][wildcards.dataset]["region"]

def get_dataset_regions_comma_separated(wildcards):
    return config["analysis_regions"]["dataset" + wildcards.number]["region"].replace(" ", ",")

def get_n_individuals(wildcards):
    return config["analysis_regions"]["simulated_dataset" + wildcards.number]["n_individuals"]

def get_n_variants(wildcards):
    return config["analysis_regions"]["simulated_dataset" + wildcards.number]["n_variants"]

def get_dataset_skipped_individuals_comma_separated(wildcards):
    return config["analysis_regions"]["svdataset" + wildcards.number]["skip_individuals"].replace(" ", ",")

rule download_sv_vcf:
    output:
        vcf="data/sv-variants.vcf.gz",
        index="data/sv-variants.vcf.gz.tbi"
    shell:
        "wget -O {output.vcf} {config[sv_vcf]} && wget -O {output.index} {config[sv_vcf]}.tbi "


rule download_vcf:
    output:
        vcf="data/variants.vcf.gz",
        index="data/variants.vcf.gz.tbi"

    shell:
        "wget -O {output.vcf} {config[thousand_genomes_vcf]} && wget -O {output.index} {config[thousand_genomes_vcf]}.tbi "


rule prepare_simulated_dataset_vcf:
    input:
        ref="data/simulated_dataset{number,\d+}/ref.fa",
        ref_index="data/simulated_dataset{number,\d+}/ref.fa.fai"
    output:
        vcf="data/simulated_dataset{number,\d+}/variants.vcf",
        vcf_gz="data/simulated_dataset{number,\d+}/variants.vcf.gz",
        individual_vcf="data/simulated_dataset{number}/truth_seed1.vcf",
        individual_vcf_gz="data/simulated_dataset{number}/truth_seed1.vcf.gz",
    params:
        n_individuals=get_n_individuals,
        n_variants=get_n_variants
    conda: "envs/prepare_data.yml"
    shell:
        "graph_read_simulator simulate_population_vcf -r {input.ref} -n {params.n_variants} -i {params.n_individuals} -o {output.vcf} -I {output.individual_vcf} && "
        "bgzip -c -f {output.vcf} > {output.vcf_gz} && tabix -f -p vcf {output.vcf_gz} && "
        "bgzip -c -f {output.individual_vcf} > {output.individual_vcf_gz} && tabix -f -p vcf {output.individual_vcf_gz} "


rule prepare_dataset_vcf:
    input:
        vcf="data/variants.vcf.gz"
    output:
        "data/dataset{number,\d+}/variants.vcf.gz"
    params:
        regions=get_dataset_regions_comma_separated
    conda: "envs/prepare_data.yml"
    shell:
        #"bcftools view -O z --regions {config[analysis_regions][{dataset}]} variants.vcf.gz > {output} && tabix -p vcf {output} "
        "bcftools view -O z --regions {params.regions} {input.vcf} > {output} && tabix -f -p vcf {output} "


rule prepare_svdataset_vcf:
    input:
        vcf="data/sv-variants.vcf.gz"
    output:
        variants="data/svdataset{number,\d+}/variants.vcf.gz",
        tmp="data/svdataset{number,\d+}/variants.vcf.gz.tmp"
    params:
        regions=get_dataset_regions_comma_separated,
        skip_individuals=get_dataset_skipped_individuals_comma_separated
    conda: "envs/prepare_data.yml"
    shell:
        #"bcftools view -O z --regions {config[analysis_regions][{dataset}]} variants.vcf.gz > {output} && tabix -p vcf {output} "
        "zcat {input.vcf} | python3 scripts/filter_uncertain_sv.py | bgzip -c -f | "
        "bcftools annotate --rename-chrs resources/chromosome-mappings.txt -O z - > {output.tmp} && tabix -p vcf -f {output.tmp} && "
        "bcftools view -f PASS --samples ^{params.skip_individuals} -O z --regions {params.regions} {output.tmp}  | "
        "bcftools +fill-tags -O z - -- -t AF "
        " > {output.variants} && tabix -f -p vcf {output.variants} "

rule prepare_dataset_reference:
    input:
        "data/hg38_chr1-Y.fa"
    output:
        full_reference="data/{dataset}/ref.fa",
        index="data/{dataset}/ref.fa.fai"
    params:
        regions=get_dataset_regions
    conda: "envs/prepare_data.yml"
    shell:
        "samtools faidx {input} {params.regions} | python3 scripts/format_fasta_header.py > {output.full_reference}  && " 
        "samtools faidx {output.full_reference}"

rule make_flat_reference:
    input: "data/{d}/ref.fa"
    output:
        fasta="data/{d}/ref_flat.fa",
        index="data/{d}/ref_flat.fa.fai"
    #shell: "grep -v '^>' {input} > {input}.no-headers && echo -e \">ref\n$(cat {input}.no-headers)\" > {output.fasta} && samtools faidx {output.fasta}"
    #shell: "grep -v '^>' {input} > {output.fasta}.tmp && echo '>ref' > {output.fasta} && cat {output.fasta}.tmp >> {output.fasta} && samtools faidx {output.fasta}"
    conda: "envs/prepare_data.yml"
    shell:
        r"""python3 scripts/make_flat_reference.py {input} | sed -e 's/.\{{80\}}/&\n/g' > {output.fasta} && samtools faidx {output.fasta}"""


rule remove_genotype_info:
    input: "{sample}.vcf.gz"
    output: "{sample}_no_genotypes.vcf"
    shell: "zcat {input} | cut -f 1-9 - > {output}"

rule get_all_sample_names_from_vcf_chinese_subpopulation:
    input:
        "resources/sample_names_chinese.txt"
    output:
        sample_names="data/{dataset}/sample_names_chinese.txt",
        sample_names_random="data/{dataset}/sample_names_random_order_chinese.txt"
    conda: "envs/prepare_data.yml"
    shell:
        "cp {input} {output.sample_names} && "
        "python scripts/shuffle_lines.py {output.sample_names} {config[random_seed]} > {output.sample_names_random} "



rule get_all_sample_names_from_vcf:
    input:
        "data/{dataset}/variants.vcf.gz"
    output:
        sample_names="data/{dataset}/sample_names_all.txt",
        sample_names_random="data/{dataset}/sample_names_random_order_all.txt"
    conda: "envs/prepare_data.yml"
    shell:
        "bcftools query -l {input} > {output.sample_names} && "
        "python scripts/shuffle_lines.py {output.sample_names} {config[random_seed]} > {output.sample_names_random} "


rule create_vcf_with_subsample_of_individuals:
    input:
        vcf="data/{dataset}/variants.vcf.gz",
        sample_names_random="data/{dataset}/sample_names_random_order_{subpopulation}.txt"
    output:
        subsamples="data/{dataset}/sample_names_random_order_{n_individuals,\d+}{subpopulation,[a-z]+}.txt",
        vcf="data/{dataset}/variants_{n_individuals,\d+}{subpopulation,[a-z]+}.vcf.gz",
        vcfindex="data/{dataset}/variants_{n_individuals,\d+}{subpopulation,[a-z]+}.vcf.gz.tbi"
    conda: "envs/prepare_data.yml"
    shell:
        "head -n {wildcards.n_individuals} {input.sample_names_random} > {output.subsamples} && "
        "bcftools view -O z -S {output.subsamples} {input.vcf} > {output.vcf} && tabix -f -p vcf {output.vcf}"


rule uncompress_subsampled_vcf:
    input:
        vcf="data/{dataset}/variants_{n_individuals,\d+}{subpopulation}.vcf.gz"
    output:
        vcf="data/{dataset}/variants_{n_individuals,\d+}{subpopulation}.vcf"
    shell:
        "zcat {input} > {output}"


rule bwa_index_reference_genome:
    input:
        "data/{dataset}/ref.fa"
    output:
        "data/{dataset}/ref.fa.bwt"
    conda: "envs/bwa.yml"
    shell:
        "bwa index {input}"


rule convert_fa_to_fq:
    input: "{reads}.fa"
    output: "{reads}.fq"
    conda: "envs/prepare_data.yml"
    shell: "scripts/convert_fa_to_fq.sh {input} > {output}"


rule make_dict_file:
    input: "data/{dataset}/ref.fa"
    output: "data/{dataset}/ref.dict"
    conda: "envs/picard.yml"
    shell: "picard CreateSequenceDictionary -R {input} -O {output}"
    