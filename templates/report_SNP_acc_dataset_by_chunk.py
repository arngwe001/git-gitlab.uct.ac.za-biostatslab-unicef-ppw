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

def acc_by_maf(inSNP_acc, outSNP_acc):
    """
    :return:
    """
    datas = {}

    outSNP_acc_out = open(outSNP_acc, 'w')
    outSNP_acc_out.writelines('\\t'.join(['Population', 'MAF<1%', 'MAF 1-5%', 'MAF>=5%'])+'\\n')
    info_datas = open(inSNP_acc).readlines()
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
            # idx_conc = data.index("concord_type0")
            idx_conc = data.index("r2_type0")
        else:
            maf = float(data[idx_exp_freq_a1])
            acc = float(data[idx_conc])
            datas[dataset]['total'] += 1
            if maf >= 0.5:
                maf = 1-maf
            if maf == 0 :
                datas[dataset]['mono'].append(acc)
            if maf > 0 and maf < 0.005:
                datas[dataset]['rare_extreme'].append(acc)
            if maf < 0.01 and maf > 0.005:
                datas[dataset]['rare'].append(acc)
            if maf > 0.05:
                datas[dataset]['common'].append(acc)
            if maf <= 0.05 and maf >= 0.01:
                datas[dataset]['moderate'].append(acc)
    for dataset in sorted(datas):
        tot = datas[dataset]['total']
        try:
            mono_ = sum(datas[dataset]['mono'])/float(len(datas[dataset]['mono']))
        except:
            mono_ = 0
        try:
            rare_extreme_ = sum(datas[dataset]['rare_extreme'])/float(len(datas[dataset]['rare_extreme']))
        except:
            rare_extreme_ = 0
        try:
            rare_ = sum(datas[dataset]['rare'])/float(len(datas[dataset]['rare']))
        except:
            rare_ = 0
        try:
            common_ = sum(datas[dataset]['common'])/float(len(datas[dataset]['common']))
        except:
            common_ = 0
        try:
            moderate_ = sum(datas[dataset]['moderate'])/float(len(datas[dataset]['moderate']))
        except:
            moderate_ = 0
        outSNP_acc_out.writelines('\\t'.join([dataset, str(format(rare_, '0,.2f')), str(format(moderate_, '0,.2f')), str(format(common_, '0,.2f'))])+'\\n')
    outSNP_acc_out.close()

args.outSNP_acc = "${outSNP_acc}"
args.inSNP_acc = "${inSNP_acc}"
args.infoCutoff = float("${infoCutoff}")
if args.inSNP_acc and args.outSNP_acc:
    acc_by_maf(args.inSNP_acc, args.outSNP_acc)
