

def get_dataset_regions_comma_separated(wildcards):
    return config["analysis_regions"][wildcards.dataset]["region"].replace(" ", ",")

# GLIMPSE has no working conda or anything
# to make installation backwards compatible, download static binaries from release
rule install_glimpse:
    output:
        "glimpse/GLIMPSE_chunk_static",
        "glimpse/GLIMPSE_phase_static",
        "glimpse/GLIMPSE_ligate_static",
    shell:
        "wget -O glimpse/GLIMPSE_chunk_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static && "
        "wget -O glimpse/GLIMPSE_phase_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_phase_static && "
        "wget -O glimpse/GLIMPSE_ligate_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_ligate_static && "
        "chmod a+x glimpse/GLIMPSE_*_static "


rule make_glimpse_chunks:
    input:
        glimpse_command="glimpse/GLIMPSE_chunk_static",
        variants="data/dataset{number}/variants.vcf.gz"
    output:
        "data/dataset{number}/glimpse_chunks.txt"
    params:
        regions=get_dataset_regions_comma_separated,
    shell:
        """
        chromosomes='{params.regions}'
        for chromosome in $(echo $chromosomes | tr "," "\n")
            do
            {input.glimpse_command} --input {input.variants} --region $chromosome  --window-size 2000000 --buffer-size 200000 --output data/dataset{wildcards.number}/glimpse_chunk.$chromosome.txt &
        done
        wait
        cat data/dataset{wildcards.number}/glimpse_chunk.*.txt > {output}
        """


rule run_glimpse:
    input:
        vcf="data/{dataset}/kageNoPriorsN{n_individuals}all_{experiment}.vcf.gz",
        ref_vcf="data/{dataset}/variants_{n_individuals}all.vcf.gz",
        glimpse_command="glimpse/GLIMPSE_phase_static",
        glimpse_command_ligate="glimpse/GLIMPSE_ligate_static",
        chunks="data/{dataset}/glimpse_chunks.txt"
    output:
        vcf="data/{dataset}/kageWithGlimpseN{n_individuals,\d+}all_{experiment}.vcf.gz"
    benchmark: "data/{dataset}/run_glimpse_{n_individuals,\d+}all_{experiment}.tsv"
    conda: "envs/bcftools.yml"
    threads: config["n_threads"]
    params:
        regions=get_dataset_regions_comma_separated,
    shell:
        """
        tabix -p vcf -f {input.vcf}
        rm -f data/{wildcards.dataset}/{input.vcf}-GLIMPSE-*.bcf
        cat {input.chunks} | parallel -j 40 --line-buffer "scripts/run_glimpse.sh {{}} {input.vcf} {input.ref_vcf}"
       
        # merge all result files 
        chromosomes='{params.regions}'
        for chromosome in $(echo $chromosomes | tr "," "\n")
            do
            LST=data/{wildcards.dataset}/glimpse_list$chromosome.tmp.txt
            ls {input.vcf}-GLIMPSE-$chromosome.*.bcf | python3 scripts/sort_glimpse_list.py > $LST
            {input.glimpse_command_ligate} --input $LST --output data/{wildcards.dataset}/glimpse_tmp_$chromosome.vcf.gz
            tabix -p vcf -f data/{wildcards.dataset}/glimpse_tmp_$chromosome.vcf.gz
        done
        wait
        bcftools concat -O z data/{wildcards.dataset}/glimpse_tmp_*.vcf.gz > {output.vcf}
        tabix -p vcf -f {output.vcf}
        
        """
    