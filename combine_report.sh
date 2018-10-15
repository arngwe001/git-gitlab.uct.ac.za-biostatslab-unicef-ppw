#!/bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:ppn=10
#PBS -j oe
#PBS -o /researchdata/fhgfs/mamana/AFRICA_CHIP/qc_two_samples/LOG/qc_two_samples.out
#PBS -M mbymam001@myuct.ac.za
#PBS -m ae

HOMEDIR="${HOME}/imputation/combine_report"
OUTDIR="/researchdata/fhgfs/mamana/combine_report"
PROJECT="gwen"

## Compile results
nextflow -log ${OUTDIR}/${PROJECT}/LOG/combine_report.log \
    run ${HOMEDIR}/combine_report.nf \
    -c ${HOMEDIR}/combine_report.config \
    -w ${OUTDIR}/work \
    -resume \
    -with-trace ${OUTDIR}/${PROJECT}/LOG/combine_report.txt \
    -with-timeline ${OUTDIR}/${PROJECT}/LOG/combine_report.html \
    -with-dag ${OUTDIR}/${PROJECT}/LOG/combine_report.dot \
    -profile pbs