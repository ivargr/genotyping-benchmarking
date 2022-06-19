


cat all_chunks.txt | parallel -j 16 --line-buffer "../../scripts/run_glimpse.sh {} bayestyper_hg002_real_reads_15x.vcf.gz.inf-replaced.vcf.gz"



