#!/bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:ppn=10
#PBS -j oe
#PBS -o /researchdata/fhgfs/arngwe001/imputation/nextflow/LOG/qc_two_samples.out
#PBS -M arngwe001001@myuct.ac.za
#PBS -m ae

HOMEDIR="${HOME}/imputation/combine_report"
OUTDIR="/researchdata/fhgfs/arngwe001/combine_report"
PROJECT="gwen"

## Compile results
nextflow -log ${OUTDIR}/${PROJECT}/LOG/combine_report.log \
    run ${HOMEDIR}/Git/combine_report.nf \
    -c ${HOMEDIR}/Git/combine_report.config \
    -w ${OUTDIR}/work \
    -resume \
    -with-trace ${OUTDIR}/${PROJECT}/LOG/combine_report.txt \
    -with-timeline ${OUTDIR}/${PROJECT}/LOG/combine_report.html \
    -with-dag ${OUTDIR}/${PROJECT}/LOG/combine_report.dot \
    -profile pbs