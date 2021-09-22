rule download_reference_genome:
    output:
        "data/hg38.2bit"
    shell:
        "wget -O {output} {config[hg38_2bit_file]}"

rule convert_reference_genome_to_fasta:
    input:
        "data/hg38.2bit"
    output:
        "data/hg38.fa"
    shell:
        "twoBitToFa {input} {output}"

rule convert_reference_to_numeric:
    input:
        "data/hg38.fa"
    output:
        "data/hg38_numeric.fa"
    shell:
        "sed 's/chr//g' {input} > {output}"

rule index_fasta:
    input:
        "data/{ref}.fa"
    output:
        "data/{ref}.fa.fai"
    shell:
        "samtools faidx -o {output} {input}"

rule remove_scaffolds_from_reference:
    input:
        fa="data/hg38_numeric.fa",
        index="data/hg38_numeric.fa.fai"
    output:
        "data/hg38_chr1-Y.fa"
    shell:
        "samtools faidx {input.fa} {config[chromosomes]} > {output}"


rule make_decoy_fasta:
    output: "data/{dataset}/decoy.fasta"
    shell: "echo -n '' > {output}"