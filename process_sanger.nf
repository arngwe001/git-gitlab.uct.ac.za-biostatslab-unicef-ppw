#!/usr/bin/env nextflow

/*
 * Authors:
 *      Mamana Mbiyavanga
 *
 *  On behalf of the H3ABionet Consortium
 *  2017
 *
 *
 * Description  : Nextflow pipeline to compile results from imputation for One population, multiple reference panels, multiple tagFiles/chip
 *
*/

//---- General definitions --------------------------------------------------//

version = '1.0'

println "=================================================="


def helpMessage() {
    log.info"""
        ${version}
    """.stripIndent()
}

// Show help emssage
//if (params.help){
//    helpMessage()
//    exit 0
//}

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
            "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

println "|-- Project : $workflow.projectDir"
println "|-- Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "|-- Command line: $workflow.commandLine"
//println "|-- Datasets: ${file(params.bedFile).getName()}, ${file(params.bimFile).getName()}, ${file(params.famFile).getName()}"


// Check if chromosomes
if (params.chromosomes == '' || params.chromosomes == 'ALL'){
    chromosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')
}
else{
    chromosomes = params.chromosomes.split(',')
}

println "|-- Chromosomes used: ${chromosomes.join(', ')}"
println "|-- Population groups used: ${params.datasets.keySet().join(', ')}"

datasets = params.datasets.keySet().join('-')
chrms = chromosomes[0]+"-"+chromosomes[-1]

// Help functions

//chromosomes = ['22']

// check if files exist
datasets_all = [:]

all_ref_names = params.ref_panels.split(',')
ref_names = params.ref_panels.split(',').join('_')


imputes = [:]
infos = [:]
refs_frq = [:]
frq_pops = [:]
imputes_ld_ = []
infos_refs_dataset = [:]
imputes_refs_dataset = [:]
infos_refs_tagName = [:]
infos_refs_dataset_all = [:]
imputes_refs_dataset_all = [:]
not_found = []
imputes_ld_all = []
chunks = []
chunk_dataset = [:]
chunk_tagName = [:]
infos_dataset = [:]
imputes_dataset = [:]
imputes_datasets = []
params.datasets.each{ dat ->
    dataset = dat.key
    all_ref_names.each { ref_name ->
        info = sprintf(dat.value.info, [ref_name,ref_name])
        impute = sprintf(dat.value.impute, [ref_name, ref_name])
        if (!file(impute).exists()) {
            not_found << impute
            exit 1, "${not_found.size()} File(s) ${not_found.join(', ')} not found. Please check your config file."
        } else {
            if (!(dataset in imputes_dataset.keySet())){
                imputes_dataset[dataset] = [dataset, ref_name+'=='+impute]
            }
            else {
                imputes_dataset[dataset][1] = imputes_dataset[dataset][1] + ',' + ref_name + '==' + impute
            }
            imputes_datasets << dataset+ '==' +ref_name+ '==' +impute
        }
    }
}


//imputes_dataset_cha = Channel.from( imputes_dataset.values() ).view()
//imputes_dataset_cha.into{ imputes_dataset_cha; imputes_dataset_1 }

"""
Fill AF and MAF in vcf
"""
process fill_vcf {
    tag "fill_${dataset}_${ref_name}"
    memory { 10.GB * task.attempt }
    publishDir "${vcfFile_path}", overwrite: true, mode:'copy'
    input:
        val impute from imputes_datasets
    output:
        set dataset, ref_name, vcfFile_path, file(vcfFile_new) into fill_vcf
    script:
        data = impute.split("==")
        dataset = data[0]
        ref_name = data[1]
        imputeFile = data[2]
        vcfFile = file(imputeFile).getName()
        vcfFile_path = file(imputeFile).getParent()
        infoFile = "${file(vcfFile).baseName}.info"
        vcfFile_new = "${file(file(vcfFile).baseName).baseName}.AF.vcf.gz"
        """
        bcftools +fill-tags ${imputeFile} -Oz -o ${vcfFile_new} -- -t AC,AN,AF,MAF
        """
}


"""
Info from VCF
"""
fill_vcf.into{fill_vcf; fill_vcf_1; fill_vcf_2}
process info_from_vcf {
    tag "info_from_vcf_${dataset}_${ref_name}"
    memory { 10.GB * task.attempt }
    publishDir "${vcfFile_path}", overwrite: true, mode:'copy'
    input:
        set dataset, ref_name, vcfFile_path, file(vcfFile) from fill_vcf_1
    output:
        set dataset, ref_name, file(infoFile) into info_from_vcf
    script:
        infoFile = "${file(file(file(vcfFile).baseName).baseName).baseName}.imputed.info"
        template "info_from_vcf.py"
}



//"""
//Filter grouped by tagname
//"""
//infos_dataset_cha = Channel.from( infos_dataset.values() )
//infos_dataset_cha.into{ infos_dataset_cha; infos_dataset_1 }
//process filter_infos_dataset {
//    tag "filter_${dataset}"
//    memory { 10.GB * task.attempt }
//    publishDir "${params.output_dir}/INFOS/${dataset}", overwrite: true, mode:'copy'
//    input:
//        set dataset, infos from infos_dataset_1
//    output:
//        set val(dataset), file(well_out), file(acc_out), file("${well_out}_snp") into filter_info_dataset
//    script:
//        chrms = chromosomes[0]+"-"+chromosomes[-1]
//        comb_info = "${dataset}_${chrms}_${params.infoCutoff.replace(".", "")}.imputed_info"
//        well_out = "${comb_info}_well_imputed"
//        acc_out = "${comb_info}_accuracy"
//        """
//        python2.7 ${params.homedir}/templates/report__.py \
//            --infoFiles ${infos} \
//            --outWell_imputed ${well_out} \
//            --outSNP_acc ${acc_out} \
//            --infoCutoff ${params.infoCutoff}
//
//        """
//}
//
//"""
//Report 1: Well imputed by dataset by chunk
//"""
//filter_info_dataset.into{ filter_info_dataset; filter_info_dataset_1 }
//process report_well_imputed_dataset {
//    tag "report_wellImputed_${dataset}_${ref_names}_${chrms}"
//    memory { 2.GB * task.attempt }
//    publishDir "${params.output_dir}/REPORTS/CSV/${dataset}", overwrite: true, mode:'copy'
//    echo true
//    input:
//        set val(dataset), file(inWell_imputed), file(acc_in), file(well_SNP) from filter_info_dataset_1
//    output:
//        set val(dataset), file(outWell_imputed) into report_well_imputed_dataset
//    script:
//        chrms = chromosomes[0]+"-"+chromosomes[-1]
//        infoCutoff = params.infoCutoff
//        outWell_imputed = "${inWell_imputed.baseName}_report_well_imputed.tsv"
//        template "report_well_imputed_dataset_by_maf.py"
//}
//
//"""
//Plot performance all tags by maf
//"""
//report_well_imputed_dataset.into{ report_well_imputed_dataset; report_well_imputed_dataset_1 }
//process plot_performance_dataset{
//    tag "plot_performance_dataset_${dataset}_${ref_names}_${chrms}"
//    cpus { 2 * task.attempt }
//    memory { 2.GB * task.cpus }
//    publishDir "${params.output_dir}/REPORTS/PLOTS/${dataset}", overwrite: true, mode:'copy'
//    input:
//        set val(dataset), file(well_imputed_report) from report_well_imputed_dataset_1
//    output:
//        set val(dataset), file(performance_by_maf_plot) into plot_performance_dataset
//    script:
//        performance_by_maf_plot = "${well_imputed_report.baseName}_performance_by_maf.tiff"
//        group = "REFERENCE_PANEL"
//        template "plot_performance_by_maf.R"
//}
//
//
//
//"""
//
//"""
//infos_dataset_cha.into{ infos_dataset_cha; infos_dataset_2 }
//process plot_r2_frequency {
//    tag "plot_r2_freq_${dataset}_${ref_names}_${chrms}"
//    publishDir "${params.output_dir}/REPORTS/PLOTS/${dataset}", overwrite: true, mode:'copy'
//    label "medium"
//    input:
//        set dataset, infos from infos_dataset_2
//    output:
//        set dataset, file(plot_out) into plot_r2_freq
//    script:
//        chrms = chromosomes[0]+"-"+chromosomes[-1]
//        plot_out = "${dataset}_${ref_names}_${chrms}_r2_freq.png"
//        impute_info_cutoff = params.infoCutoff
//        template "r2_Frequency_plot.R"
//}


//"""
//Report 2: Accuracy by dataset by chunk
//"""
//filter_info_dataset.into{ filter_info_dataset; filter_info_dataset_2}
//process report_SNP_acc_dataset {
//    tag "report_SNP_acc_${dataset}_${all_ref_names}_${chrms}"
//    memory { 2.GB * task.attempt }
//    maxRetries 1
//    publishDir "${params.output_dir}/REPORTS/CSV/${dataset}", overwrite: true, mode:'copy'
//    input:
//        set val(dataset), file(well_in), file(inSNP_acc) from filter_info_dataset_2
//    output:
//        set val(dataset), file(outSNP_acc) into report_SNP_acc_dataset
//    script:
//        chrms = chromosomes[0]+"-"+chromosomes[-1]
//        infoCutoff = params.infoCutoff
//        outSNP_acc = "${inSNP_acc.baseName}_report_SNP_acc.tsv"
//        template "report_SNP_acc_dataset_by_chunk.py"
//}
//
//"""
//Plot accuracy all tags by maf
//"""
//report_SNP_acc_dataset.into{ report_SNP_acc_dataset; report_SNP_acc_dataset_1 }
//process plot_report_SNP_acc_dataset{
//    tag "plot_SNP_acc_dataset_${dataset}"
//    cpus { 2 * task.attempt }
//    memory { 2.GB * task.cpus }
//    publishDir "${params.output_dir}/REPORTS/PLOTS/${dataset}", overwrite: true, mode:'copy'
//    input:
//        set val(dataset), file(SNP_acc_report) from report_SNP_acc_dataset_1
//    output:
//        set val(dataset), file(SNP_acc_by_maf_plot) into plot_report_SNP_acc_dataset
//    script:
//        SNP_acc_by_maf_plot = "${SNP_acc_report.baseName}_SNP_acc_by_maf.tiff"
//        group = "CHIP"
//        template "plot_SNP_acc_by_maf.R"
//}
//
//
//"""
//Report: Well imputed for each tagName by maf
//"""
//filter_info_tagName.into{ filter_info_tagName; filter_info_tagName_1}
//process report_well_imputed_tagName {
//    tag "report_wellImputed_${tagName}_${datasets}_${chrms}"
//    memory { 10.GB * task.attempt }
//    publishDir "${params.output_dir}/REPORTS/CSV/${tagName}", overwrite: true, mode:'copy'
//    input:
//        set val(tagName), file(inWell_imputed), file(acc_in) from filter_info_tagName_1
//    output:
//        set val(tagName), file(outWell_imputed) into report_well_imputed_tagName
//    script:
//        chrms = chromosomes[0]+"-"+chromosomes[-1]
//        infoCutoff = params.impute_info_cutoff
//        comb_info = "${params.dataset}_${tagName}_${all_ref_names}_${chrms}_${infoCutoff.replace(".", "")}.imputed_info"
//        outWell_imputed = "${comb_info}_report_well_imputed.tsv"
//        template "report_well_imputed_dataset_by_maf.py"
//}
//
//
//"""
//Plot performance all populations for each tag by maf
//"""
////report_well_imputed_tagName.into{ report_well_imputed_tagName; report_well_imputed_tagName_1 }
////process plot_performance_tagName{
////    tag "plot_performance_tagName_${tagName}_${datasets}_${chrms}"
////    cpus { 2 * task.attempt }
////    memory { 2.GB * task.cpus }
////    publishDir "${params.output_dir}/REPORTS/PLOTS/${tagName}", overwrite: true, mode:'copy'
////    input:
////        set val(tagName), file(well_imputed_report) from report_well_imputed_tagName_1
////    output:
////        set val(tagName), file(performance_by_maf_plot) into plot_performance_tagName
////    script:
////        chrms = chromosomes[0]+"-"+chromosomes[-1]
////        performance_by_maf_plot = "${well_imputed_report.baseName}_performance_by_maf.tiff"
////        group = "CLUSTER"
////        template "plot_performance_by_maf.R"
////}
//
//
//
//"""
//Report 2.1: Accuracy for each tagName for all population by maf
//"""
//filter_info_tagName.into{ filter_info_tagName; filter_info_tagName_2}
//process report_SNP_acc_tagName {
//    tag "report_SNP_acc_${tagName}_${datasets}_${chrms}"
//    memory { 10.GB * task.attempt }
//    maxRetries 1
//    publishDir "${params.output_dir}/REPORTS/CSV/${tagName}", overwrite: true, mode:'copy'
//    input:
//        set val(tagName), file(well_in), file(inSNP_acc) from filter_info_tagName_2
//    output:
//        set val(tagName), file(outSNP_acc) into report_SNP_acc_tagName
//    script:
//        chrms = chromosomes[0]+"-"+chromosomes[-1]
//        infoCutoff = params.impute_info_cutoff
//        comb_info = "${params.dataset}_${tagName}_${all_ref_names}_${chrms}_${infoCutoff.replace(".", "")}.imputed_info"
//        outSNP_acc = "${comb_info}_report_SNP_acc.tsv"
//        template "report_SNP_acc_dataset_by_chunk.py"
//}
///
//"""
///Plot accuracy all populations by maf
//"""
////report_SNP_acc_tagName.into{ report_SNP_acc_tagName; report_SNP_acc_tagName_1 }
////process plot_report_SNP_acc_tagName{
///    tag "plot_SNP_acc_tagName_${tagName}_${datasets}_${chrms}"
////    cpus { 2 * task.attempt }
////    memory { 2.GB * task.cpus }
////    publishDir "${params.output_dir}/REPORTS/PLOTS/${tagName}", overwrite: true, mode:'copy'
////    input:
////        set val(tagName), file(SNP_acc_report) from report_SNP_acc_tagName_1
////    output:
////        set val(tagName), file(SNP_acc_by_maf_plot) into plot_report_SNP_acc_tagName
////    script:
////        SNP_acc_by_maf_plot = "${SNP_acc_report.baseName}_SNP_acc_by_maf.tiff"
////        group = "CLUSTER"
////        template "plot_SNP_acc_by_maf.R"
////}
//
//
//// """
//// Report: Coverage per chunk
//// """
//// infos_refs_cha = Channel.from(infos_refs_dataset.values())
//// process report_cov_chunk {
////    tag "report_${dataset}_${chrm}_${tags}"
////    memory { 20.GB * task.attempt }
////    publishDir "${params.output_dir}/REPORTS/CSV/${dataset}", overwrite: true, mode:'copy'
////    input:
////        set val(dataset), val(chrm), val(chunks) from infos_refs_cha
////    output:
////        set val(dataset), val(chrm), file(report_out) into report_cov_chunk
////    script:
////        tags = params.tagSNPs_files.keySet().join('-')
////        report_out = "${dataset}_${chrm}_${tags}_report.tsv"
////        """
////        python2.7 -u ${params.homedir}/templates/report__.py \
////           --dataset ${dataset} \
////           --chunks \'${chunks}\' \
////           --cov_chunk \'${report_out}\' \
////           --infoCutoff ${params.impute_info_cutoff} \
////           --tagNames \'${params.tagSNPs_files.keySet().join('==')}\' \
////           --info_format ${params.output_dir}/IMPUTE_RESULTS/{tagName}/${dataset}/${chrm}/${dataset}_{tagName}_${chrm}_{chunk_start}_{chunk_end}_imputed_info \
////           --pop_maf_format ${params.output_dir}/IMPUTE_RESULTS/{tagName}/${dataset}/${chrm}/FRQ/${dataset}_{tagName}_chr${chrm}_{chunk_start}_{chunk_end}_${params.ref_1.name}-${params.ref_2.name}.frq \
////           --ref_maf_format ${params.output_dir}/IMPUTE_RESULTS/${dataset}/${chrm}/FRQ/${dataset}_${chrm}_{chunk_start}_{chunk_end}_merged.frq
////
////        """
//// }
////
//// report_cov_chunk.into{ report_cov_chunk; report_cov_chunk_1 }
//// process plot_results{
////     tag "plot_results_${dataset}_${chrm}_${tags}"
////     cpus { 2 * task.attempt }
////     memory { 2.GB * task.cpus }
////     publishDir "${params.output_dir}/REPORTS/PLOTS/${dataset}/${chrm}", overwrite: true, mode:'copy'
////     input:
////         set val(dataset), val(chrm), file(report_csv) from report_cov_chunk_1
////     output:
////         set val(dataset), val(chrm), file(coverage_plot), file(well_imputed_plot), file(average_r2_plot), file(average_concordance_plot), file(common_snps_on_chip_plot) into plot_results
////     script:
////         tags = params.tagSNPs_files.keySet().join('-')
////         coverage_plot = "${report_csv.baseName}_coverage.tiff"
////         well_imputed_plot = "${report_csv.baseName}_well_imputed.tiff"
////         average_r2_plot = "${report_csv.baseName}_average_r2.tiff"
////         average_concordance_plot = "${report_csv.baseName}_average_concordance_plot.tiff"
////         common_snps_on_chip_plot = "${report_csv.baseName}_common_snps_on_chip_plot.tiff"
////         template "plot_coverage_per_pop.R"
//// }
//
//"""
//Combine all chunk infos into chromosome by dataset
//"""
//infos_refs_dataset_all_cha = Channel.from( infos_refs_dataset_all.values() )
//process info_Combine_to_chrom {
//    tag "infoComb_${dataset}_${tagName}_${chrm}"
//    memory { 2.GB * task.attempt }
//    publishDir "${params.output_dir}/INFOS/${dataset}/${tagName}/", overwrite: true, mode:'copy'
//    input:
//        set val(dataset), val(tagName), val(chrm), val(info_files) from infos_refs_dataset_all_cha
//    output:
//        set val(dataset), val(tagName), val(chrm), file(comb_info) into info_Combine_to_chrom
//    script:
//        comb_info = "${params.dataset}_${dataset}_${tagName}_${chrm}_${params.impute_info_cutoff.replace(".", "")}.imputed_info"
//        """
//        ## Info files
//        echo "snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0" > ${comb_info}
//        tail -q -n +2 ${info_files} >> ${comb_info}
//        """
//}
//info_Combine_to_chrom.into { info_Combine_to_chrom; info_Combine_to_chrom_1 }
//info_Combine_to_chrom_list = info_Combine_to_chrom_1.toSortedList().val
//dataset_infos_all = [:]
//info_Combine_to_chrom_list.each{ dataset, tagName, chrm, comb_info ->
//    // if (!(tagName in tag_infos.keySet())){
//    //     tag_infos[tagName] = [tagName, dataset+"=="+comb_info]
//    // }
//    // else{
//    //     tag_infos[tagName][1] = tag_infos[tagName][1] + ',' + dataset+"=="+comb_info
//    // }
//    if (!(dataset in dataset_infos_all.keySet())){
//        dataset_infos_all[dataset] = [dataset, tagName+"=="+chrm+"=="+comb_info]
//    }
//    else{
//        dataset_infos_all[dataset][1] += ',' + tagName+"=="+chrm+"=="+comb_info
//    }
//}
//
//
//"""
//Report: per chromosome
//"""
//dataset_infos_all_cha = Channel.from(dataset_infos_all.values())
//process report_performance_chrm {
//   tag "report_${dataset}_${chrms}_${tags}"
//   memory { 20.GB * task.attempt }
//   publishDir "${params.output_dir}/REPORTS/CSV/${dataset}", overwrite: true, mode:'copy'
//   input:
//       set val(dataset), val(infos) from dataset_infos_all_cha
//   output:
//       set val(dataset), file(report_out) into report_performance_chrm
//   script:
//       report_out = "${params.dataset}_${dataset}_${tags}_${chrms}_${params.impute_info_cutoff.replace(".", "")}_report.tsv"
//       """
//       python2.7 -u ${params.homedir}/templates/report__.py \
//          --dataset ${dataset} \
//          --infos \'${infos}\' \
//          --cov_chrm \'${report_out}\' \
//          --infoCutoff ${params.impute_info_cutoff}
//       """
//}
//
//
////report_performance_chrm.into{ report_performance_chrm; report_performance_chrm_1 }
////process plot_report_performance_chrm{
////    tag "plot_performance_chrm_${dataset}_${tags}"
////    cpus { 2 * task.attempt }
////    memory { 2.GB * task.cpus }
////    publishDir "${params.output_dir}/REPORTS/PLOTS/${dataset}", overwrite: true, mode:'copy'
////    input:
////        set val(dataset), file(report_csv) from report_performance_chrm_1
////    output:
////        set val(dataset), file(well_imputed_plot), file(average_r2_plot), file(average_concordance_plot) into plot_report_performance_chrm
////    script:
////        well_imputed_plot = "${report_csv.baseName}_well_imputed.tiff"
////        average_r2_plot = "${report_csv.baseName}_average_r2.tiff"
////        average_concordance_plot = "${report_csv.baseName}_average_concordance_plot.tiff"
////        template "plot_coverage_per_pop_chrm.R"
////}
//
//
//
//"""
//Combine all chunk imputed data into one by dataset
//dataset refers to refName
//"""
//process imputeCombine {
//    tag "impComb_${refName}_${tagName}_${chrm}"
//    memory { 2.GB * task.attempt }
//    publishDir "${params.impute_result}/combined", overwrite: true, mode:'copy'
//    input:
//        set val(refName), val(tagName), val(chrm), val(imputed_files) from imputes_cha
//    output:
//        set val(refName), val(tagName), val(chrm), file(comb_impute) into imputeCombine
//    script:
//        comb_impute = "${params.project}_${refName}_${tagName}_chr${chrm}.imputed.gz"
//        """
//        zcat ${imputed_files} | bgzip -c > ${comb_impute}
//        """
//}
//
//
//dataset_sample_cha = Channel.from(datasets_all.values())
//process dataset_sample {
//    tag "generate_dataset_sample_${dataset}"
//    memory { 2.GB * task.attempt }
//    publishDir "${params.output_dir}/sample_files", overwrite: true, mode:'copy'
//    input:
//        set val(dataset), val(dataset_vcfs) from dataset_sample_cha
//    output:
//        set val(dataset), file("${dataset}.sample") into dataset_sample
//    script:
//        dataset_vcfs_ = dataset_vcfs.split(' ')
//        dataset_sample = "${file(file(dataset_vcfs_[0]).baseName).baseName}.sample"
//        """
//        bcftools \
//            query --list-samples \
//            ${dataset_vcfs_[0]} \
//            > ${dataset_sample}
//        python2.7 ${params.homedir}/templates/update_sample.py \
//            --inSample ${dataset_sample} \
//            --outSample ${dataset}.sample
//        """
//}
//
//imputeCombine.into { imputeCombine; imputeCombine_1 }
//dataset_sample_list = dataset_sample.toSortedList().val
//filter_info_dataset_snp.into{ filter_info_dataset_snp; filter_info_dataset_snp_1 }
//filter_info_dataset_snp_1_list = filter_info_dataset_snp_1.toSortedList().val
//def transform_imputeCombine = { refName, dataset, chrm, comb_impute  ->
//    imputeCombine_sample = []
//    dataset_sample_list.each{ dataset_, sample_file ->
//        if(dataset == dataset_){
//            filter_info_dataset_snp_1_list.each{ refName_, well_snp_file ->
//                if(refName == refName_) {
//                    imputeCombine_sample << [refName, dataset, chrm, comb_impute, sample_file, well_snp_file]
//                }
//            }
//        }
//    }
//    return imputeCombine_sample
//}
//
//imputeCombine_sample_cha = imputeCombine_1
//        .flatMap{ it -> transform_imputeCombine(it) }
//
//process imputeToVCF {
//    tag "toVCF_${refName}_${tagName}_${chrm}"
//    memory { 2.GB * task.attempt }
//    publishDir "${params.impute_result}/plink", overwrite: true, mode:'copy'
//    input:
//        set refName, tagName, val(chrm), file(imputed_gz), file(dataset_sample), file(well_snp_file) from imputeCombine_sample_cha
//    output:
//        set refName, tagName, val(chrm), file("${imputed_gz.baseName}.vcf.gz") into imputeToVCF
//    script:
//        """
//        gunzip -c ${imputed_gz} > ${imputed_gz.baseName}
//        plink2 \
//            --gen ${imputed_gz.baseName} \
//            --sample ${dataset_sample} \
//            --oxford-single-chr ${chrm} \
//            --hard-call-threshold 0.1 \
//            --allow-no-sex \
//            --extract ${well_snp_file} \
//            --recode vcf-iid --out ${imputed_gz.baseName} || true
//        bgzip ${imputed_gz.baseName}.vcf
//        rm -f ${imputed_gz.baseName}
//        """
//}
//
//imputeToVCF.into{imputeToVCF; imputeToVCF_1}
//imputeToVCF_1_list = imputeToVCF_1.toSortedList().val
//imputeCombine_chrm_all = [:]
//imputeToVCF_1_list.each{ refName, dataset, chrm, vcfFile  ->
//    if(!(refName in imputeCombine_chrm_all)){
//        imputeCombine_chrm_all[refName] = [refName, dataset, vcfFile]
//    }
//    else{
//        imputeCombine_chrm_all[refName][2] += ' '+vcfFile
//    }
//}
//
//imputeCombine_chrm_all_cha = Channel.from(imputeCombine_chrm_all.values())
//
//process concat_dataset_vcf {
//    tag "concat_${dataset}_${refName}_chr${chrms}"
//    memory { 50.GB * task.attempt }
//    publishDir "${params.output_dir}/POP_DATA/VCF", overwrite: true, mode: 'symlink'
//    input:
//        set val(refName), val(dataset), val(dataset_vcfs) from imputeCombine_chrm_all_cha
//    output:
//        set val(refName), val(dataset), file(vcf_out) into concat_dataset_vcf
//    script:
//        vcf_out = "${dataset}_${refName}_chr${chrms}.vcf.gz"
//        """
//        bcftools concat \
//            ${dataset_vcfs} \
//            -Oz -o ${dataset}.tmp.vcf.gz
//        ## Recalculate AC, AN, AF
//        bcftools +fill-tags  ${dataset}.tmp.vcf.gz -Oz -o  ${dataset}.tmp1.vcf.gz
//        bcftools sort ${dataset}.tmp1.vcf.gz -Oz -o ${vcf_out}
//        rm ${dataset}.tmp*.vcf.gz
//        """
//}
//        tags = params.tagSNPs_files.keySet().join('-')
//        chrms = chromosomes[0]+"-"+chromosomes[-1]
//        infoCutoff = params.impute_info_cutoff
//        outWell_imputed = "${inWell_imputed.baseName}_report_well_imputed.tsv"
//        template "report_well_imputed_dataset_by_maf.py"
//
//"""
//Plot performance all tags by maf
//"""
////report_well_imputed_dataset.into{ report_well_imputed_dataset; report_well_imputed_dataset_1 }
////process plot_performance_dataset{
////    tag "plot_performance_dataset_${dataset}_${tags}_${chrms}"
////    cpus { 2 * task.attempt }
////    memory { 2.GB * task.cpus }
////    publishDir "${params.output_dir}/REPORTS/PLOTS/${dataset}", overwrite: true, mode:'copy'
////    input:
////        set val(dataset), file(well_imputed_report) from report_well_imputed_dataset_1
////    output:
////        set val(dataset), file(performance_by_maf_plot) into plot_performance_dataset
////    script:
////        performance_by_maf_plot = "${well_imputed_report.baseName}_performance_by_maf.tiff"
////        group = "CHIP"
////        template "plot_performance_by_maf.R"
////}
//
//
//"""
//Repor 2: Accuracy by dataset by chunk
//"""
//filter_info_dataset.into{ filter_info_dataset; filter_info_dataset_2}
//process report_SNP_acc_dataset {
//    tag "report_SNP_acc_${dataset}_${tags}_${chrms}"
//    memory { 2.GB * task.attempt }
//    maxRetries 1
//        set val(refName), val(dataset), val(dataset_vcfs) from imputeCombine_chrm_all_cha
//        set val(refName), val(dataset), file(vcf_out) into concat_dataset_vcf
//        vcf_out = "${dataset}_${refName}_chr${chrms}.vcf.gz"
//        """
//        bcftools concat \
//            ${dataset_vcfs} \
//            -Oz -o ${dataset}.tmp.vcf.gz
//        ## Recalculate AC, AN, AF
//        bcftools +fill-tags  ${dataset}.tmp.vcf.gz -Oz -o  ${dataset}.tmp1.vcf.gz
//        bcftools sort ${dataset}.tmp1.vcf.gz -Oz -o ${vcf_out}
//        rm ${dataset}.tmp*.vcf.gz
//        """

////workflow.onComplete {
////    if ( workflow.success == true ){
////        def subject = 'My pipeline execution'
////        def recipient = 'mbymam001@myuct.ac.za'
////
////        ['mail', '-s', subject, recipient].execute() << """
////
////        Pipeline execution summary
////        ---------------------------
////        Commandline: ${workflow.commandLine}
////        Configuration file(s): ${workflow.configFiles}
////        Completed at: ${workflow.complete}
////        Duration    : ${workflow.duration}
////        Success     : ${workflow.success}
////        workDir     : ${workflow.workDir}
////        exit status : ${workflow.exitStatus}
////        Error report: ${workflow.errorReport ?: '-'}
////        """
////  }
////}
