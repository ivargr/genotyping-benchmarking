# Modified version of Pangenie run-from-callset snakemake pipeline (https://bitbucket.org/jana_ebler/pangenie/raw/d340042099a3af7683729327235a3cd3c8bb8a71/pipelines/run-from-callset/Snakefile)

# assign IDs to all alleles
rule add_ids:
	input:
		"data/{dataset}/variants_{n_individuals}individuals.vcf.gz"
	output:
		"data/{dataset}/variants_{n_individuals,\d+}individuals_with_id.vcf"
	shell:
		'zcat {input} | python3 scripts/pangenie_add_ids.py > {output}'


# merge variants into a pangenome graph
rule merge_haplotypes:
	input:
		vcf='data/{dataset}/variants_{n_individuals}individuals_with_id.vcf',
		reference="data/{dataset}/ref.fa"
	output:
		"data/{dataset}/variants_{n_individuals}individuals_multiallelic.vcf"
	shell:
		"""
		python3 scripts/pangenie_merge_vcfs.py merge -vcf {input.vcf} -r {input.reference} -ploidy 2 > {output}
		"""
