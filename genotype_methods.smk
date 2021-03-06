


rule run_kmc:
    input:
        "data/{dataset}/{reads}.fa"
    output:
        pre="data/{dataset}/{reads}.kmc.out.kmc_pre",
        suf="data/{dataset}/{reads}.kmc.out.kmc_suf"
    benchmark:
        "data/{dataset}/benchmarks/malva_kmc_{reads}.tsv"
    threads:
        config["n_threads"]
    resources:
        mem_gb=20
    conda:
        "envs/kmc.yml"
    shell:
        "rm -rf kmc_tmp && "
        "mkdir -p kmc_tmp && "
        "kmc -ci1 -m8 -k43 -t{config[n_threads]} -fm {input} data/{wildcards.dataset}/{wildcards.reads}.kmc.out kmc_tmp"


rule run_malva:
    input:
        ref="data/{dataset}/ref.fa",
        variants="data/{dataset}/variants.vcf.gz",
        kmer_pre="data/{dataset}/{reads}.kmc.out.kmc_pre",
        kmer_suf="data/{dataset}/{reads}.kmc.out.kmc_suf"
    output:
        genotypes="data/{dataset}/malva_{reads}.vcf.gz",
        benchmark="data/{dataset}/benchmarks/malva_{reads}.tsv"
    #benchmark:
    #    "data/{dataset}/benchmarks/malva_{reads}.tsv"
    conda:
        "envs/malva.yml"
    resources:
        mem_gb=100
    threads: 4
    shell:
        "/usr/bin/time -v malva-geno -k 35 -r 43 -b 16 {input.ref} {input.variants} data/{wildcards.dataset}/{wildcards.reads}.kmc.out 1> {output.genotypes}.tmp 2> {output.benchmark} && "
        # convert nan quality to 0
        """awk '{{ $6 = ($6 == "nan" ? 0 : $6) }} 1' OFS="\\t" {output.genotypes}.tmp | bgzip -f -c > {output.genotypes}"""
        #"malva-geno call -k 35 -r 41 -b 16 {input.ref} {input.variants} {wildcards.reads}.kmc.out > {output.genotypes}"

rule run_pangenie:
    input:
        ref = "data/{dataset}/ref.fa",
        reads="data/{dataset}/{reads}.fa",
        variants = "data/{dataset}/variants_{n_individuals}all_multiallelic.vcf",
    output:
        genotypes="data/{dataset}/pangenieN{n_individuals,\d+}_{reads}.vcf",
        genotypes_gz="data/{dataset}/pangenieN{n_individuals,\d+}_{reads}.vcf.gz",
        benchmark="data/{dataset}/benchmarks/pangenieN{n_individuals}_{reads}.tsv"
    threads:
        config["n_threads"]
    resources:
        mem_gb=150
    #benchmark:
        #"data/{dataset}/benchmarks/pangenie_{reads}.{n_individuals}individuals.tsv"
    #    "data/{dataset}/benchmarks/pangenieN{n_individuals}_{reads}.tsv"
    conda: "envs/bwa.yml"
    shell:
        "/usr/bin/time -v {config[pangenie_path]} -i {input.reads} -r {input.ref} -v {input.variants} -j 32 -k 31 " 
        " -t {config[n_threads]} -g -o {output.genotypes} 2> {output.benchmark} "
        "&& mv {output.genotypes}_genotyping.vcf {output.genotypes} && "
        "bgzip -f -c {output.genotypes} > {output.genotypes_gz}  && tabix -f -p vcf {output.genotypes_gz}"


# hack to run pangenie with 32 individuals if none are specified
rule pangenie_wrapper:
    input:
        genotypes = "data/{dataset}/pangenieN32_{experiment}.vcf.gz"
    output:
        genotypes = "data/{dataset}/pangenie_{experiment,[a-z0-9_]+}.vcf.gz"
    conda: "envs/bwa.yml"
    shell:
        "cp {input} {output} && tabix -f -p vcf {output}"


rule map_reads_to_reference:
    input:
        ref="data/{dataset}/ref.fa",
        index="data/{dataset}/ref.fa.bwt",
        reads="data/{dataset}/{sample_id,\w+}_{read_info}.fq"
    output:
        bam="data/{dataset}/{sample_id,\w+}_{read_info}.mapped.bam",
        sorted_bam="data/{dataset}/{sample_id,\w+}_{read_info}.mapped.sorted.bam",
        sorted_bam_index="data/{dataset}/{sample_id,\w+}_{read_info}.mapped.sorted.bam.bai"
    threads:
        config["n_threads"]
    resources:
        mem_gb=80
    benchmark:
        "data/{dataset}/benchmarks/bwamem_{sample_id}_{read_info}.tsv"
    conda: "envs/bwa.yml"
    shell:
        "bwa mem -R '@RG\\tID:{wildcards.sample_id}\\tSM:{wildcards.sample_id}' -t {config[n_threads]} {input.ref} {input.reads} | samtools view -b > {output.bam} && "
        "sambamba sort {output.bam} && sambamba index -p {output.sorted_bam}"


rule run_graphtyper:
    input:
        sorted_bam="data/{dataset}/{reads}.mapped.sorted.bam",
        fai="data/{dataset}/ref.fa.fai",
        ref= "data/{dataset}/ref.fa",
        variants = "data/{dataset}/variants.vcf.gz",
    output:
        genotypes="data/{dataset}/graphtyper_{reads}.vcf.gz",
        benchmark="data/{dataset}/benchmarks/graphtyper_{reads}.tsv"
    threads:
        config["n_threads"]
    #benchmark:
    #    "data/{dataset}/benchmarks/graphtyper_{reads}.tsv"
    conda: "envs/graphtyper.yml"
    shell:
        "/usr/bin/time -v ./scripts/run_graphtyper.sh {wildcards.reads} {input.ref} {input.sorted_bam} {input.variants} {output.genotypes} {config[n_threads]} 2> {output.benchmark}"
        #"rm -rf graphtyper_results_{wildcards.reads} && "
        #"python3 scripts/run_graphtyper_in_parallel.py 2500000 {input.fai} "
        #"| parallel --line-buffer -j {config[n_threads]} 'graphtyper genotype {input.ref} --output=graphtyper_results_{wildcards.reads} --sam={input.sorted_bam} --region={{}} --threads=1 --vcf {input.variants}' && "
        #"find graphtyper_results_{wildcards.reads}/* -name '*.vcf.gz' | sort > vcf_file_list && "
        #"bcftools concat --naive --file-list vcf_file_list -Oz | bcftools sort -Oz > {output.genotypes}"


rule run_kmc_bayestyper:
    input:
        "data/{dataset}/{reads}.fa"
    output:
        pre="data/{dataset}/{reads}.kmers_bayestyper.kmc_pre",
        suf="data/{dataset}/{reads}.kmers_bayestyper.kmc_suf"
    threads:
        config["n_threads"]
    resources:
        mem_gb=20
    benchmark:
        "data/{dataset}/benchmarks/bayestyper_kmc_{reads}.tsv"
    conda:
        "envs/kmc.yml"
    shell:
        "mkdir -p kmc_tmp_bayestyper && "
        "kmc -t{config[n_threads]} -k55 -ci1 -fa {input} data/{wildcards.dataset}/{wildcards.reads}.kmers_bayestyper kmc_tmp_bayestyper"


rule make_multiallelic_variants_for_bayestyper:
    input:
        vcf="data/{dataset}/variants_no_genotypes.vcf",
        ref= "data/{dataset}/ref.fa",
    output:
        "data/{dataset}/variants_no_genotypes_multiallelic.vcf"
    conda: "envs/bcftools.yml"
    shell:
        "bcftools norm -m +any -f {input.ref} {input.vcf} > {output}"


rule make_samples_tsv_for_bayestyper:
    output:
        "data/{dataset}/samples_{individual,\w+}_{reads_type}_{coverage,\d+}x.tsv"
    shell:
        "echo '{wildcards.individual}\tM\tdata/{wildcards.dataset}/{wildcards.individual}_{wildcards.reads_type}_{wildcards.coverage}x.kmers_bayestyper' > {output}"


rule make_bloomfilter_for_bayestyper:
    input:
        pre = "data/{dataset}/{reads}.kmers_bayestyper.kmc_pre",
        suf = "data/{dataset}/{reads}.kmers_bayestyper.kmc_suf"
    output:
        "data/{dataset}/{reads}.kmers_bayestyper.bloomData",
        "data/{dataset}/{reads}.kmers_bayestyper.bloomMeta"
    threads:
        config["n_threads"]
    benchmark:
        "data/{dataset}/benchmarks/bayestyper_bloomfilter_{reads}.tsv"
    conda: "envs/bayestyper.yml"
    shell:
        "bayesTyperTools makeBloom -k data/{wildcards.dataset}/{wildcards.reads}.kmers_bayestyper -p {config[n_threads]}"



rule run_bayestyper:
    input:
        ref="data/{dataset}/ref.fa",
        decoy="data/{dataset}/decoy.fasta",
        bloomdata="data/{dataset}/{reads}.kmers_bayestyper.bloomData",
        bloommeta="data/{dataset}/{reads}.kmers_bayestyper.bloomMeta",
        variants="data/{dataset}/variants_no_genotypes_multiallelic.vcf",
        samples="data/{dataset}/samples_{reads}.tsv",
    output:
        #units=dynamic("data/{dataset}/tmp_bayestyper_data_{reads}/bayestyper_unit_{unit_id}/variant_clusters.bin")
        #units="data/{dataset}/tmp_bayestyper_data_{reads}/bayestyper_unit_1/variant_clusters.bin"
        genotypes="data/{dataset,\w+}/bayestyper_{reads}.vcf.gz"
    benchmark:
        "data/{dataset}/benchmarks/bayestyper_{reads}.tsv"
    params:
        out_prefix ="bayestyper/bayestyper",
    threads: config["n_threads"]
    resources:
        mem_gb=50
    conda: "envs/bayestyper.yml"
    shell:
        "mkdir -p data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads} && "
        "rm -rf data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/* && "
        "bayesTyper cluster -v {input.variants} -s {input.samples} -g {input.ref} -d {input.decoy} -p {config[n_threads]} -o data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper && "
        "for dir in data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper_unit_*; do\n "
        "    bayesTyper genotype -v $dir/variant_clusters.bin -c data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper_cluster_data -s {input.samples} -g {input.ref} -d {input.decoy} -o $dir/bayestyper -z -p {config[n_threads]}; "
        "done && "
        # hack since bcftools seem to not bgzip, but only gzip:
        "ls data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper_unit_*/*.vcf.gz | xargs -P 16 -n 1 gunzip && "
        "ls data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper_unit_*/*.vcf | xargs -P 16 -n 1 bgzip && "
        "ls data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper_unit_*/*.vcf.gz | xargs -P 16 -n 1 tabix -f -p vcf && "
        
        "bcftools concat -a -O z -o {output.genotypes} data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper_unit_*/bayestyper.vcf.gz"


rule run_gatk:
    input:
        sorted_bam="data/{dataset}/{reads}.mapped.sorted.bam",
        fai="data/{dataset}/ref.fa.fai",
        ref="data/{dataset}/ref.fa",
        dict="data/{dataset}/ref.dict",
        variants="data/{dataset}/variants.vcf.gz",
    output:
        genotypes="data/{dataset}/gatk_{reads}.vcf.gz"
    threads:
        config["n_threads"]
    benchmark:
        "data/{dataset}/benchmarks/gatk_{reads}.tsv"
    resources:
        mem_gb=50
    conda: "envs/gatk.yml"
    shell:
        "rm -f data/{wildcards.dataset}/{wildcards.reads}_tmp_gatk_variants_*.vcf && "
        "rm -f data/{wildcards.dataset}/{wildcards.reads}_tmp_gatk_variants_*.vcf.gz && "
        "rm -f data/{wildcards.dataset}/{wildcards.reads}_tmp_gatk_variants_*.vcf.gz.tbi && "
        "python3 scripts/partition_genome.py 10000000 {input.fai} | "
        "parallel --line-buffer -j {config[n_threads]} gatk HaplotypeCaller "
        "--reference {input.ref} --input {input.sorted_bam} --native-pair-hmm-threads 1 --output data/{wildcards.dataset}/{wildcards.reads}_tmp_gatk_variants_{{}}.vcf --intervals {{}} "
        "--minimum-mapping-quality 20 "
        #"--alleles {input.variants} "
        "&& "
        "ls data/{wildcards.dataset}/{wildcards.reads}_tmp_gatk_variants_*.vcf | xargs -P 16 -n 1 bgzip -f && "
        "ls data/{wildcards.dataset}/{wildcards.reads}_tmp_gatk_variants_*.vcf.gz | xargs -P 16 -n 1 tabix -f -p vcf && "
        "python3 scripts/partition_genome.py 10000000 data/{wildcards.dataset}/ref.fa.fai 'data/{wildcards.dataset}/{wildcards.reads}_tmp_gatk_variants_---.vcf.gz' > gatk_file_list && "
        "bcftools concat -a -O z --file-list gatk_file_list | bcftools sort -O z > {output.genotypes}"
