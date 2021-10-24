# Modified version of Pangenie run-from-callset snakemake pipeline (https://bitbucket.org/jana_ebler/pangenie/raw/d340042099a3af7683729327235a3cd3c8bb8a71/pipelines/run-from-callset/Snakefile)


## parameters
vcf = config['vcf']
reference = config['reference']
scripts = 'scripts'
# allow max 20 % missing alleles at a variant position
frac_missing = 1.0
outdir = config['outdir']
pangenie = config['pangenie']
samples = config['reads'].keys()

###########################################################
##  Convert VCF into pangenome graph by merging variants ##
##  that are overlapping into multi-allelic positions.   ##
###########################################################

# check that VCF is correct and remove positions with more than frac_missing missing alleles
rule prepare_vcf:
	input:

	output:
		"data/{dataset}/variants_"
		temp('{outdir}/input-vcf/input-missing-removed.vcf')
	log:
		'{outdir}/input-vcf/prepare-vcf.log'
	shell:
		"cat {input} | python3 {scripts}/prepare-vcf.py --missing {frac_missing} 2> {log} 1> {output}"


# assign IDs to all alleles
rule add_ids:
	input:
		'{outdir}/input-vcf/input-missing-removed.vcf'
	output:
		'{outdir}/input-vcf/callset.vcf'
	log:
		'{outdir}/input-vcf/callset.log'
	shell:
		'cat {input} | python3 {scripts}/add-ids.py 2> {log} 1> {output}'


# create biallelic VCF with one record per ALT allele
rule normalize:
	input:
		'{outdir}/input-vcf/callset.vcf'
	output:
		'{outdir}/input-vcf/callset-biallelic.vcf'
	log:
		'{outdir}/input-vcf/callset-biallelic.log'
	shell:
		'bcftools norm -m- {input} 2> {log} 1> {output}'


# merge variants into a pangenome graph
rule merge_haplotypes:
	input:
		vcf = '{outdir}/input-vcf/callset-biallelic.vcf',
		reference = reference
	output:
		'{outdir}/pangenome/pangenome.vcf'
	log:
		 '{outdir}/pangenome/pangenome.log'
	resources:
		mem_total_mb=900000,
		runtime_hrs=30,
		runtime_min=1
	shell:
		"""
		python3 {scripts}/merge_vcfs.py merge -vcf {input.vcf} -r {input.reference} -ploidy 2  2> {log} 1> {output}
		"""


##########################################################
##                  Run PanGenie                        ##
##########################################################

rule pangenie:
	input:
		vcf='{outdir}/pangenome/pangenome.vcf',
		reference = reference,
		reads = lambda wildcards: config['reads'][wildcards.sample]
	output:
		'{outdir}/pangenie/{sample}_graph_genotyping.vcf'
	threads:
		24
	resources:
		mem_total_mb=200000,
		runtime_hrs=4,
		runtime_min=59
	log:
		'{outdir}/pangenie/pangenie-{sample}.log'
	benchmark:
		'{outdir}/benchmarks/pangenie-{sample}.txt'
	params:
		prefix = "{outdir}/pangenie/{sample}_graph"
	shell:
		"{pangenie} -i {input.reads} -v {input.vcf} -r {input.reference} -o {params.prefix} -j {threads} -t {threads} -g &> {log}"



##########################################################
##     Convert VCF back to original representation      ##
##########################################################

# represent genotypes in a bi-allelic VCF with one record per ALT-allele
rule convert_back_biallelic_representation:
	input:
		pangenie='{outdir}/pangenie/{sample}_graph_genotyping.vcf',
		panel='{outdir}/input-vcf/callset.vcf'
	output:
		'{outdir}/pangenie/{sample}_genotypes-biallelic.vcf'
	resources:
		mem_total_mb=20000
	benchmark:
		'{outdir}/benchmarks/{sample}-convert-to-biallelic.txt'
	shell:
		'cat {input.pangenie} | python3 {scripts}/convert-to-biallelic.py {input.panel} > {output}'

rule compress_vcf:
	input:
		"{filename}.vcf"
	output:
		"{filename}.vcf.gz"
	shell:
		"""
		bgzip -c {input} > {output}
		tabix -p vcf {output}
		"""


# represent genotypes in the same way as in the input callset
rule convert_back_original_representation:
	input:
		pangenie='{outdir}/pangenie/{sample}_genotypes-biallelic.vcf.gz'
	output:
		'{outdir}/genotypes/{sample}-genotypes.vcf'
	benchmark:
		'{outdir}/benchmarks/{sample}-convert-to-original.txt'
	resources:
		mem_total_mb=10000
	conda:
		"env/merging.yml"
	shell:
		'bcftools norm -m+ {input.pangenie} > {output}'