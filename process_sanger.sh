#!/bin/bash
#PBS -q UCTlong
#PBS -l nodes=1:ppn=10
#PBS -j oe
#PBS -o /researchdata/fhgfs/arngwe001/imputation/nextflow/LOG/qc_two_samples.out
#PBS -M arngwe001@myuct.ac.za
#PBS -m ae

HOMEDIR="${HOME}/imputation/gwen"
OUTDIR="/researchdata/fhgfs/mamana/process_sanger"
PROJECT="gwen"
#ls ${HOMEDIR}/process_sanger.nf
## Compile results
cd ${OUTDIR}
nextflow -log ${OUTDIR}/${PROJECT}/LOG/process_sanger.log \
    run ${HOMEDIR}/process_sanger.nf \
    -c ${HOMEDIR}/combine_report.config \
    -w ${OUTDIR}/work \
#    -resume \
#    -with-trace ${OUTDIR}/${PROJECT}/LOG/process_sanger.txt \
#    -with-timeline ${OUTDIR}/${PROJECT}/LOG/process_sanger.html \
#    -with-dag ${OUTDIR}/${PROJECT}/LOG/process_sanger.dot \
#    -profile pbs


