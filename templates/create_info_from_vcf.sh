#!/bin/bash
FILE=$1
echo $FILE
echo -e "SNP\tREF(0)\tALT(1)\tALT_Frq\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose0\tDose1" > ${FILE}.info
echo annotate
bcftools +fill-tags ${FILE}.gz -Oz -o ${FILE}_FREQ.vcf.gz -- -t AF,MAF
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' ${FILE}_FREQ.vcf.gz -Oz -o ${FILE}_SNPID.vcf.gz
echo query
bcftools query -f '%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/INFO\t1\t-\t-\t-\t-\t-\t-\n' ${FILE}_SNPID.vcf.gz >> ${FILE}.info
bgzip -f ${FILE}.info
