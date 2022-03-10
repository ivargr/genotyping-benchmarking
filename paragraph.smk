


rule make_paragraph_samples_file:
    input:
        sorted_bam="data/{dataset}/{sample_id}_{read_info}_{coverage}x.mapped.sorted.bam",
    output:
        sample_file="data/{dataset}/paragraph_samples_{sample_id,[a-zA-Z0-9]+}_{read_info}_{coverage,\d+}x.txt"
    shell:
        """
        echo 'id\tpath\tdepth\tread length\n{wildcards.sample_id}\t{input.sorted_bam}\t{wildcards.coverage}\t{config[read_length]}' > {output.sample_file}
        """
        
        
rule run_paragraph:
    input:
        variants="data/{dataset}/variants.vcf.gz",
        sample_file="data/{dataset}/paragraph_samples_{sample_id,\w+}_{read_info}.txt",
        ref="data/{dataset}/ref.fa"
    output:
        genotypes="data/{dataset}/paragraph_{sample_id,\w+}_{read_info}/genotypes.vcf.gz"
    shell:
        """
        python3 /home/ivar/dev/paragraph/bin/multigrmpy.py -i {input.variants} -m {input.sample_file} -r {input.ref} -o data/{wildcards.dataset}/paragraph_{wildcards.sample_id}_{wildcards.read_info}/ -t {config[n_threads]}
        """


rule postprocess_paragraph:
    input:
        genotypes = "data/{dataset}/paragraph_{sample_id}_{reads}/genotypes.vcf.gz"
    output:
        genotypes = "data/{dataset}/paragraph_{sample_id,[a-zA-Z0-9]+}_{reads}.vcf.gz",
        index = "data/{dataset}/paragraph_{sample_id,[a-zA-Z0-9]+}_{reads}.vcf.gz.tbi"
    shell:
        "bcftools view -s {wildcards.sample_id} -O z {input.genotypes} > {output.genotypes} && "
        "tabix -p vcf -f {output.genotypes}"
