#!/usr/bin/env python2.7

import gzip,time
import argparse,sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("--vcfFile", help="")
parser.add_argument("--infoFile", help="")
args = parser.parse_args()

def vcf_info(vcf, infoOUT):
    info = open(infoOUT, 'w')
    info.writelines("\\t".join(["SNP", "REF(0)",  "ALT(1)",  "ALT_Frq", "MAF",     "AvgCall", "Rsq",     "Genotyped",       "LooRsq",  "EmpR",    "EmpRsq",  "Dose0",   "Dose1"])+"\\n")
    for line in gzip.open(vcf):
        if not line.startswith('##'):
            data = line.strip().split('\t')
            if "#CHROM" in line:
                CHROM_idx = data.index("#CHROM")
                POS_idx = data.index("POS")
                ID_idx = data.index("ID")
                REF_idx = data.index("REF")
                ALT_idx = data.index("ALT")
                INFO_idx = data.index("INFO")
            else:
                if "TYPED" in line:
                    GEN = "Genotyped"
                else:
                    GEN = "Imputed"
                    INFO = data[INFO_idx].split(";")
                    INFOS = {}
                    for inf in INFO:
                        dat = inf.split('=')
                        try:
                            INFOS[dat[0]] = dat[1]
                        except:
                            print dat
                info.writelines("\\t".join([data[CHROM_idx][1:]+'_'+data[POS_idx], data[REF_idx], data[ALT_idx], INFOS['AF'], INFOS['MAF'], '0.99', INFOS['INFO'], GEN, '-', '-', '-', '-', '-'])+"\\n")
                # time.sleep(5)
    info.close()

args.vcfFile="${vcfFile}"
args.infoFile="${infoFile}"
vcf_info(args.vcfFile, args.infoFile)