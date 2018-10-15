#!/bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:ppn=10
#PBS -j oe
#PBS -o /researchdata/fhgfs/mamana/AFRICA_CHIP/qc_two_samples/LOG/qc_two_samples.out
#PBS -M mbymam001@myuct.ac.za
#PBS -m ae

. $HOME/.bashrc


HOMEDIR="${HOME}/chipdesign/imputation_chip_design"
OUTDIR="/researchdata/fhgfs/mamana/AFRICA_CHIP"
PROJECT="qc_three_samples"

## load Python virtual environment
source activate imputation


## Nextflow script here
mkdir -p "${OUTDIR}/chipdesign_nf_1/"
cd "${OUTDIR}/chipdesign_nf_1/"
#nextflow -log ${OUTDIR}/${PROJECT}/LOG/imputation_impute2_2refs.log \
#    run ${HOMEDIR}/imputation_impute2_1.nf \
#    -c ${HOMEDIR}/HPC/imputation_impute2_2refs.config \
#    -w ${OUTDIR}/work \
#    -resume \
#    -with-trace ${OUTDIR}/${PROJECT}/LOG/imputation_impute2_2refs.txt \
#    -with-timeline ${OUTDIR}/${PROJECT}/LOG/imputation_impute2_2refs.html \
#    -with-dag ${OUTDIR}/${PROJECT}/LOG/imputation_impute2_2refs.dot \
#    -profile pbs

## Compile results
nextflow -log ${OUTDIR}/${PROJECT}/LOG/combine_report_2.log \
    run ${HOMEDIR}/combine_report_2.nf \
    -c ${HOMEDIR}/HPC/imputation_impute2_2refs.config \
    -w ${OUTDIR}/work \
    -resume \
    -with-trace ${OUTDIR}/${PROJECT}/LOG/combine_report_2.txt \
    -with-timeline ${OUTDIR}/${PROJECT}/LOG/combine_report_2.html \
    -with-dag ${OUTDIR}/${PROJECT}/LOG/combine_report_2.dot \
    -profile pbs