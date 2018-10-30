#!/bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:ppn=10
#PBS -j oe
#PBS -o /researchdata/fhgfs/arngwe001/imputation/nextflow/LOG/qc_two_samples.out
#PBS -M arngwe001@myuct.ac.za
#PBS -m ae

HOMEDIR="${HOME}/imputation/Git"
OUTDIR="/researchdata/fhgfs/arngwe001/combine_report"
PROJECT="gwen"
#ls ${HOMEDIR}/combine_report.nf
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