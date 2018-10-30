#!/bin/bash
for FIL in ASW_AGR_chr22.imputed EUR_AGR_chr22.imputed  GAMBIA_AGR_chr22.imputed  KENYA_AGR_chr22.imputed  ZULU_AGR_chr22.imputed
do
    ~/imputation/Git/templates/create_info_from_vcf.sh ${FIL}
done