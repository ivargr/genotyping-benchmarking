


cat chunks.*.txt > all_chunks.txt

#/usr/bin/time cat chunks.1.txt | parallel -j 16 --line-buffer "../../scripts/run_glimpse.sh {} kageNoHelperModelN250all_hg002_simulated_reads_15x.vcf.gz"
#cat all_chunks.txt | parallel -j 16 --line-buffer "../../scripts/run_glimpse.sh {} gatk_hg002_real_reads_30x.vcf.gz.inf-replaced.vcf.gz"
cat all_chunks.txt | parallel -j 16 --line-buffer "../../scripts/run_glimpse.sh {} bayestyper_hg002_real_reads_15x.vcf.gz.inf-replaced.vcf.gz"



