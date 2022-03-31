


rule make_plink_map:
    input:
        "data/{dataset}/variants.vcf.gz"
    output:
        "data/{dataset}/variants.map"
    shell:
        "plink --vcf {input} --recode --out data/{wildcards.dataset}/variants"


rule refine_genotypes_with_beagle:
    input:
        variants="data/{dataset}/gatk_{reads}.vcf.gz",
        ref_variants="data/{dataset}/variants.vcf.gz",
        plink_map="data/{dataset}/variants.map"
    output:
        refined_variants="data/{dataset}/gatkRefined_{reads}.vcf.gz"
    conda:
        "envs/beagle4.yml"
    threads:
        config["n_threads"]
    shell:
        "beagle gl={input.variants} ref={input.ref_variants} out={output.refined_variants} map={input.plink_map} impute=false nthreads={config[n_threads]}"

