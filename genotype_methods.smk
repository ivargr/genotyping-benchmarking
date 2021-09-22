


rule run_kmc:
    input:
        "data/{dataset}/{reads}.fa"
    output:
        pre="data/{dataset}/{reads}.kmc.out.kmc_pre",
        suf="data/{dataset}/{reads}.kmc.out.kmc_suf"
    shell:
        "kmc -ci0 -m8 -k41 -t{config[n_threads]} -fm {input} data/{wildcards.dataset}/{wildcards.reads}.kmc.out ."

rule run_malva:
    input:
        ref="data/{dataset}/ref.fa",
        variants="data/{dataset}/variants_no_overlaps.vcf.gz",
        kmer_pre="data/{dataset}/{reads}.kmc.out.kmc_pre",
        kmer_suf="data/{dataset}/{reads}.kmc.out.kmc_suf"
    output:
        genotypes="data/{dataset}/malva_{reads}.vcf"
    shell:
        "malva-geno -k 35 -r 41 -b 16 {input.ref} {input.variants} data/{wildcards.dataset}/{wildcards.reads}.kmc.out > {output.genotypes}.tmp && "
        # convert nan quality to 0
        """awk '{{ $6 = ($6 == "nan" ? 0 : $6) }} 1' OFS="\\t" {output.genotypes}.tmp > {output.genotypes}"""
        #"malva-geno call -k 35 -r 41 -b 16 {input.ref} {input.variants} {wildcards.reads}.kmc.out > {output.genotypes}"

"""
rule run_malva_genotyping:
    input:
        kmers_pre="data/{dataset}/{reads}.kmc.out.kmc_pre",
        kmers_suf="data/{dataset}/{reads}.kmc.out.kmc_suf",
        variants="data/{dataset}/variants_no_overlaps_no_genotypes.vcf",
        ref="data/{dataset}/ref.fa"

    output:
        "data/{dataset}/malva_genotypes_{reads}.vcf"
"""

rule run_pangenie:
    input:
        ref = "data/{dataset}/ref.fa",
        reads="data/{dataset}/{reads}.fa",
        variants = "data/{dataset}/variants_no_overlaps_{n_individuals}individuals.vcf",
    output:
        genotypes="data/{dataset}/pangenie_{reads}.{n_individuals,\d+}individuals.vcf",
        #genotypes_gz="data/{dataset}/pangenie_{reads}.{n_individuals,\d+}individuals.vcf.gz"
    threads:
        config["n_threads"]
    shell:
        "PanGenie -i {input.reads} -r {input.ref} -v {input.variants} -j 32 -k 31 -t {config[n_threads]} -g -o {output.genotypes} && mv {output.genotypes}_genotyping.vcf {output.genotypes} "
        #"bgzip -f -c {output.genotypes} > {output.genotypes_gz} && tabix -p vcf {output.genotypes_gz}"


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
    shell:
        "bwa mem -R '@RG\\tID:{wildcards.sample_id}\\tSM:{wildcards.sample_id}' -t {config[n_threads]} {input.ref} {input.reads} | samtools view -b > {output.bam} && "
        "sambamba sort {output.bam} && sambamba index -p {output.sorted_bam}"


rule run_graphtyper:
    input:
        sorted_bam="data/{dataset}/{reads}.mapped.sorted.bam",
        fai="data/{dataset}/ref.fa.fai",
        ref= "data/{dataset}/ref.fa",
        variants = "data/{dataset}/variants_no_overlaps.vcf.gz",
    output:
        genotypes="data/{dataset}/graphtyper_{reads}.vcf.gz"
    threads:
        config["n_threads"]
    shell:
        "rm -rf results && "
        "python3 scripts/run_graphtyper_in_parallel.py 2500000 {input.fai} "
        "| parallel --line-buffer -j {config[n_threads]} 'graphtyper genotype {input.ref} --sam={input.sorted_bam} --region={{}} --threads=1 --vcf {input.variants}' && "
        "bcftools concat --naive results/*/*.vcf.gz > {output.genotypes}"


rule run_kmc_bayestyper:
    input:
        "data/{dataset}/{reads}.fa"
    output:
        pre="data/{dataset}/{reads}.kmers_bayestyper.kmc_pre",
        suf="data/{dataset}/{reads}.kmers_bayestyper.kmc_suf"
    shell:
        "kmc -t{config[n_threads]} -k55 -ci1 -fa {input} data/{wildcards.dataset}/{wildcards.reads}.kmers_bayestyper ."


rule make_multiallelic_variants_for_bayestyper:
    input:
        vcf="data/{dataset}/variants_no_overlaps_no_genotypes.vcf",
        ref= "data/{dataset}/ref.fa",
    output:
        "data/{dataset}/variants_no_overlaps_no_genotypes_multiallelic.vcf"
    shell:
        "bcftools norm -m +any -f {input.ref} {input.vcf} > {output}"


rule make_samples_tsv_for_bayestyper:
    output:
        "data/{dataset}/samples_{individual}_simulated_reads_{coverage,\d+}x.tsv"
    shell:
        "echo '{wildcards.individual}\tM\tdata/{wildcards.dataset}/{wildcards.individual}_simulated_reads_{wildcards.coverage}x.kmers_bayestyper' > {output}"


rule make_bloomfilter_for_bayestyper:
    input:
        pre = "data/{dataset}/{reads}.kmers_bayestyper.kmc_pre",
        suf = "data/{dataset}/{reads}.kmers_bayestyper.kmc_suf"
    output:
        "data/{dataset}/{reads}.kmers_bayestyper.bloomData",
        "data/{dataset}/{reads}.kmers_bayestyper.bloomMeta"
    shell:
        "bayesTyperTools makeBloom -k data/{wildcards.dataset}/{wildcards.reads}.kmers_bayestyper -p {config[n_threads]}"



rule run_bayestyper:
    input:
        ref="data/{dataset}/ref.fa",
        decoy="data/{dataset}/decoy.fasta",
        bloomdata="data/{dataset}/{reads}.kmers_bayestyper.bloomData",
        bloommeta="data/{dataset}/{reads}.kmers_bayestyper.bloomMeta",
        variants="data/{dataset}/variants_no_overlaps_no_genotypes_multiallelic.vcf",
        samples="data/{dataset}/samples_{reads}.tsv",
    output:
        #units=dynamic("data/{dataset}/tmp_bayestyper_data_{reads}/bayestyper_unit_{unit_id}/variant_clusters.bin")
        #units="data/{dataset}/tmp_bayestyper_data_{reads}/bayestyper_unit_1/variant_clusters.bin"
        genotypes="data/{dataset,\w+}/bayestyper_{reads}.vcf.gz"
    params:
        out_prefix ="bayestyper/bayestyper",
    threads: config["n_threads"]
    shell:
        "mkdir -p data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads} && "
        "rm -rf data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/* && "
        "bayesTyper cluster -v {input.variants} -s {input.samples} -g {input.ref} -d {input.decoy} -p {config[n_threads]} -o data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper && "
        "for dir in data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper_unit_*; do\n "
        "    bayesTyper genotype -v $dir/variant_clusters.bin -c data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper_cluster_data -s {input.samples} -g {input.ref} -d {input.decoy} -o $dir/bayestyper -z -p {config[n_threads]}; "
        "done && "
        "bcftools concat -O z -o {output.genotypes} data/{wildcards.dataset}/tmp_bayestyper_data_{wildcards.reads}/bayestyper_unit_*/bayestyper.vcf.gz"


rule run_gatk:
    input:
        sorted_bam="data/{dataset}/{reads}.mapped.sorted.bam",
        fai="data/{dataset}/ref.fa.fai",
        ref="data/{dataset}/ref.fa",
        dict="data/{dataset}/ref.dict",
        variants="data/{dataset}/variants_no_overlaps.vcf.gz",
    output:
        genotypes="data/{dataset}/gatk_{reads}.vcf.gz"
    threads:
        config["n_threads"]
    shell:
        "rm -f data/{wildcards.dataset}/{wildcards.reads}_tmp_gatk_variants_*.vcf && "
        "python3 scripts/partition_genome.py 1000000 {input.fai} | "
        "parallel --line-buffer -j {config[n_threads]} /home/ivar/dev/gatk-4.1.9.0/gatk HaplotypeCaller "
        "--reference {input.ref} --input {input.sorted_bam} --output data/{wildcards.dataset}/{wildcards.reads}_tmp_gatk_variants_{{}}.vcf --intervals {{}} "
        "--minimum-mapping-quality 20 --alleles {input.variants} && "
        "bcftools concat -O z data/{wildcards.dataset}/{wildcards.reads}_tmp_gatk_variants_*.vcf > {output.genotypes}"
        #"cat gatk_genotypes.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > gatk_genotypes_sorted.vcf
