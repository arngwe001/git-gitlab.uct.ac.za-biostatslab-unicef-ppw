#!/usr/bin/env python2.7

import argparse,sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("--infoFiles", help="")
parser.add_argument("--outWell_imputed", help="")
parser.add_argument("--outSNP_acc", help="")
parser.add_argument("--infoCutoff", help="")
parser.add_argument("--inWell_imputed", help="")
parser.add_argument("--inSNP_acc", help="")
parser.add_argument("--report_acc", help="")
parser.add_argument("--ldFiles", help="")
parser.add_argument("--report_ld", help="")
parser.add_argument("--cov_chunk", help="Imputation coverage report by chunks for all tag lists for a population")
parser.add_argument("--cov_chrm", help="Imputation coverage report by chrm for all tag lists for a population")
parser.add_argument("--infos", help="")
parser.add_argument("--chunk_files", help="")
parser.add_argument("--maf_file", help="")
parser.add_argument("--dataset", help="")
parser.add_argument("--info_format", help="")
parser.add_argument("--pop_maf_format", help="")
parser.add_argument("--ref_maf_format", help="")
parser.add_argument("--tagNames", help="list of tagNames (chip list), sep by comma")
parser.add_argument("--chrms", help="list of chromosomes, sep by comma")
parser.add_argument("--chunks", help="list of chunks, sep by comma")

args = parser.parse_args()

def well_imputed_by_maf(inWell_imputed, outWell_imputed):
    """

    :return:
    """
    datas = {}
    outWell_imputed_out = open(outWell_imputed, 'w')
    outWell_imputed_out_1 = open(outWell_imputed+"_summary.tsv", 'w')
    outWell_imputed_out.writelines('\t'.join(['Population', 'MAF<1%', 'MAF 1-5%', 'MAF>=5%'])+'\\n')
    outWell_imputed_out_1.writelines('\t'.join(['Population', 'MAF<1%', 'MAF 1-5%', 'MAF>=5%', 'Total'])+'\\n')
    info_datas = open(inWell_imputed).readlines()
    for line in info_datas:
        data = line.strip().split()
        dataset = data[0]
        if dataset not in datas:
            datas[dataset] = {}
            datas[dataset]['moderate'] = []
            datas[dataset]['rare'] = []
            datas[dataset]['rare_extreme'] = []
            datas[dataset]['mono'] = []
            datas[dataset]['common'] = []
            datas[dataset]['total'] = 0
        if "snp_id" in line and "info" in line:
            idx_exp_freq_a1 = data.index('exp_freq_a1')
        else:
            maf = float(data[idx_exp_freq_a1])
            datas[dataset]['total'] += 1
            if maf > 0.5:
                maf = 1-maf
            # if maf == 0 :
            #     datas[dataset]['mono'].append(maf)
            # if maf > 0 and maf < 0.005:
            #     datas[dataset]['rare_extreme'].append(maf)
            if maf > 0.005:
                datas[dataset]['rare'].append(maf)
            if maf > 0.05:
                datas[dataset]['common'].append(maf)
            if maf <= 0.05 and maf >= 0.01:
                datas[dataset]['moderate'].append(maf)
    for dataset in sorted(datas):
        tot = datas[dataset]['total']
        # mono = len(datas[dataset]['mono'])
        # rare_extreme = len(datas[dataset]['rare_extreme'])
        rare = len(datas[dataset]['rare'])
        moderate = len(datas[dataset]['moderate'])
        common = len(datas[dataset]['common'])
        # if rare_extreme >= 1000000: rare_extreme_ = str(format(rare_extreme/1000000., '0,.1f'))+'M ('
        # else: rare_extreme_ = str(format(rare_extreme, '0,'))+' ('
        # if mono >= 1000000: mono_ = str(format(mono/1000000., '0,.1f'))+'M ('
        # else: mono_ = str(format(mono, '0,'))+' ('
        if rare >= 1000000: rare_ = str(format(rare/1000000., '0,.1f'))+'M ('
        else: rare_ = str(format(rare, '0,'))+' ('
        if moderate >= 1000000: moderate_ = str(format(moderate/1000000., '0,.1f'))+'M ('
        else: moderate_ = str(format(moderate, '0,'))+' ('
        if common >= 1000000: common_ = str(format(common/1000000., '0,.1f'))+'M ('
        else: common_ = str(format(common, '0,'))+' ('
        outWell_imputed_out.writelines('\t'.join([dataset, \
            str(len(datas[dataset]['rare'])), \
            str(len(datas[dataset]['moderate'])), \
            str(len(datas[dataset]['common']))])+'\\n')
        outWell_imputed_out_1.writelines('\t'.join([dataset, \
            rare_+str(rare * 100/tot)+'%)', \
            moderate_+str(moderate * 100/tot)+'%)', \
            common_+str(common * 100/tot)+'%)', \
            str(format(tot/1000000., '0,.1f'))+'M']) +'\\n')
    outWell_imputed_out.close()
    outWell_imputed_out_1.close()


args.inWell_imputed = "${inWell_imputed}"
args.outWell_imputed = "${outWell_imputed}"
if args.inWell_imputed and args.outWell_imputed:
    well_imputed_by_maf(args.inWell_imputed, args.outWell_imputed)
