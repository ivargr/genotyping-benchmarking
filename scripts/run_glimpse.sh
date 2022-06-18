glimpse_params=$1
ID=$(echo $glimpse_params| cut -d" " -f1)
chromosome=$(echo $glimpse_params| cut -d" " -f2)
IRG=$(echo $glimpse_params| cut -d" " -f3)
ORG=$(echo $glimpse_params | cut -d" " -f4)


REF="variants.vcf.gz"
VCF=$2
MAP="../../resources/genetic_maps/chr$chromosome.b38.gmap.gz"

OUT=$chromosome.${ID}.bcf

echo $IRG

GLIMPSE_phase_static --input-GL --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
bcftools index -f ${OUT}
