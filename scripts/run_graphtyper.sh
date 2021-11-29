reads_name=$1
ref=$2
bam=$3
variants=$4
output=$5
n_threads=$6


rm -rf graphtyper_results_$reads_name
python3 scripts/run_graphtyper_in_parallel.py 2500000 ${ref}.fai | \
parallel --line-buffer -j $n_threads "graphtyper genotype $ref --output=graphtyper_results_$reads_name --sam=$bam --region={} --threads=1 --vcf $variants"
find graphtyper_results_$reads_name/* -name '*.vcf.gz' | sort > vcf_file_list
bcftools concat --naive --file-list vcf_file_list -Oz | bcftools sort -Oz > $output