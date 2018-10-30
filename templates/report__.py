#!/usr/bin/env python2.7

import argparse,sys
import time
import gzip

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

def filter_info(infoFiles, infoCutoff, outWell_imputed, outSNP_acc):
    """
    Return:
        well_imputed: certainy >= 1
        SNP_concordance: concord_type0 != -1
    """
    well_imputed = {}
    well_imputed_common = {}
    SNP_concordance = {}
    count = 0
    infoFiles = infoFiles.split(',')
    datas = {}
    header = []
    outWell_imputed_out = open(outWell_imputed, 'w')
    outWell_imputed_snp_out = open(outWell_imputed+"_snp", 'w')
    outSNP_accuracy_out = open(outSNP_acc, 'w')
    print infoFiles
    for infoFile in infoFiles:
        infoFile = infoFile.strip().split('==')
        dataset = infoFile[0]
        info = infoFile[1]
        well_imputed[dataset] = []
        well_imputed_common[dataset] = []
        SNP_concordance[dataset] = []
        # outWell_imputed_out_dataset = open(dataset+'_'+outWell_imputed, 'w')
        for line in gzip.open(info):
            data = line.strip().split()
            if "SNP" in line and "Rsq" in line:
                if len(header) == 0:
                    header = data
                    info_idx = header.index("Rsq")
                    conc_idx = header.index("LooRsq")
                    certainty_idx = header.index("AvgCall")
                    maf_idx = header.index("MAF")
                    #pos_idx = header.index("position")
                    a1_idx = header.index('REF(0)')
                    a0_idx = header.index('ALT(1)')
                    outWell_imputed_out.writelines(' '.join([dataset]+data)+'\n')
                    # outWell_imputed_out_dataset.writelines(' '.join([dataset]+data)+'\n')
                    outWell_imputed_snp_out.writelines(data[1]+'\n')
                    outSNP_accuracy_out.writelines(' '.join([dataset]+data)+'\n')
            else:
                info = data[info_idx]
                conc = data[conc_idx]
                #pos = data[pos_idx]
                a0 = data[a0_idx]
                a1 = data[a1_idx]
                try:
                    exp_maf = float(data[maf_idx])
                except:
                    exp_maf = 0
                if exp_maf == 0:
                    maf = 0
                elif exp_maf >= 0.5:
                    maf = 1 - exp_maf
                else:
                    maf = exp_maf
                certainty = data[certainty_idx]
                if info != '-' and (float(info) >= float(infoCutoff) and float(certainty) == 1):
                    outWell_imputed_out.writelines(' '.join([dataset]+data)+'\n')
                    if maf >= 0.05:
                        well_imputed_common[dataset].append(data)
                    well_imputed[dataset].append(data)
                    outWell_imputed_snp_out.writelines(data[1]+'\n')
                if conc == '-':
                    outSNP_accuracy_out.writelines(' '.join([dataset]+data)+'\n')
                count += 1
    outWell_imputed_out.close()
    outWell_imputed_snp_out.close()
    outSNP_accuracy_out.close()
    return well_imputed[dataset],well_imputed_common[dataset]

def well_imputed_by_maf(inWell_imputed, outWell_imputed):
    """

    :return:
    """
    datas = {}
    outWell_imputed_out = open(outWell_imputed, 'w')
    outWell_imputed_out.writelines('\t'.join(['Population', 'MAF<1%', 'MAF 1-5%', 'MAF>=5%'])+'\n')
    info_datas = open(inWell_imputed).readlines()
    for line in info_datas:
        data = line.strip().split()
        dataset = data[0]
        if dataset not in datas:
            datas[dataset] = {}
            datas[dataset]['1_5'] = []
            datas[dataset]['1'] = []
            datas[dataset]['5'] = []
            datas[dataset]['total'] = 0
        if "snp_id" in line and "info" in line:
            idx_exp_freq_a1 = data.index('exp_freq_a1')
        else:
            maf = float(data[idx_exp_freq_a1])
            datas[dataset]['total'] += 1
            if maf >= 0.5:
                maf = 1-maf
            if maf >= 0.01:
                datas[dataset]['1'].append(maf)
                if maf >= 0.05:
                    datas[dataset]['5'].append(maf)
                if maf <= 0.05 and maf >= 0.01:
                    datas[dataset]['1_5'].append(maf)
    for dataset in sorted(datas):
        tot = datas[dataset]['total']
        # TODO Check this !!!!
        outWell_imputed_out.writelines('\t'.join([dataset, str(format(len(datas[dataset]['5'])/1000000., '0,.3f')), str(format(len(datas[dataset]['1'])/1000000., '0,.3f')), str(format(len(datas[dataset]['1_5'])/1000000., '0,.3f'))])+'\n')
    outWell_imputed_out.close()

def acc_by_maf(inSNP_acc, outSNP_acc):
    """
    :return:
    """
    datas = {}

    outSNP_acc_out = open(outSNP_acc, 'w')
    outSNP_acc_out.writelines('\t'.join(['Population', 'MAF<1%', 'MAF 1-5%', 'MAF>=5%'])+'\n')
    info_datas = open(inSNP_acc).readlines()
    for line in info_datas:
        data = line.strip().split()
        dataset = data[0]
        if dataset not in datas:
            datas[dataset] = {}
            datas[dataset]['1_5'] = []
            datas[dataset]['1'] = []
            datas[dataset]['5'] = []
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
            if maf >= 0.01:
                datas[dataset]['1'].append(acc)
                if maf >= 0.05:
                    datas[dataset]['5'].append(acc)
                if maf <= 0.05 and maf >= 0.01:
                    datas[dataset]['1_5'].append(acc)
    for dataset in sorted(datas):
        tot = datas[dataset]['total']
        print tot
        outSNP_acc_out.writelines('\t'.join([dataset, str(format(sum(datas[dataset]['5'])/float(len(datas[dataset]['5'])), '0,.3f')), str(format(sum(datas[dataset]['1'])/float(len(datas[dataset]['1'])), '0,.3f')), str(format(sum(datas[dataset]['1_5'])/float(len(datas[dataset]['1_5'])), '0,.3f'))])+'\n')
    outSNP_acc_out.close()


def ld_by_maf(ldFiles, report_ld, inWell_imputed, infoCutoff):
    """

    :return:
    """
    maf_data = {}
    print 'Reading', inWell_imputed
    # inWell_imputed_data = open(inWell_imputed).readlines()
    for line in open(inWell_imputed):
        data = line.split(' ')
        dataset = data[0]
        if "snp_id" in line and "info" in line:
            idx_exp_freq_a1 = data.index('exp_freq_a1')
            rs_id = data.index('rs_id')
        if dataset not in maf_data:
            maf_data[dataset] = {}
        maf_data[dataset][data[rs_id]] = data[idx_exp_freq_a1]
    datas = {}
    ld_data = {}
    report_ld_out = open(report_ld, 'w')
    report_ld_out.writelines('\t'.join(['Population', 'MAF<1%', 'MAF 1-5%', 'MAF>5%', 'TOTAL'])+'\n')
    ldFiles = ldFiles.split(',')
    not_ = 0
    in_ = 0
    header = []
    infoCutoff = float(infoCutoff)
    for ldFile in ldFiles:
        ldFile = ldFile.strip().split('==')
        dataset =ldFile[0]
        if dataset not in datas:
            datas[dataset] = {}
            datas[dataset]['1_5'] = 0
            datas[dataset]['1'] = 0
            datas[dataset]['5'] = 0
            datas[dataset]['ALL'] = set()
            datas[dataset]['total'] = len(maf_data[dataset])
        ld = ldFile[1]
        print 'Reading', ld
        ld_data[dataset] = []
        for line in open(ld):
            data = line.strip().split()
            if "SNP_A" in line and "SNP_B" in line:
                if len(header) == 0:
                    header = data
                    snpA_idx = header.index("SNP_A")
                    snpB_idx = header.index("SNP_B")
                    r2_idx = header.index("R2")
            else:
                if float(data[r2_idx]) >= infoCutoff:
                    snpA = data[snpA_idx]
                    snpB = data[snpB_idx]
                    for snp in [snpA, snpB]:
                        if snp not in datas[dataset]['ALL']:
                            try:
                                datas[dataset]['ALL'].add(snp)
                                maf = float(maf_data[dataset][snp])
                                if maf >= 0.5:
                                    maf = 1-maf
                                if maf >= 0.01:
                                    datas[dataset]['1'] += 1
                                    if maf >= 0.05:
                                        datas[dataset]['5'] += 1
                                    if maf <= 0.05 and maf >= 0.01:
                                        datas[dataset]['1_5'] += 1
                                # in_ += 1
                            except:
                                continue
                                # not_ += 1
    for dataset in sorted(datas):
        tot = datas[dataset]['total']
        # tot = len(datas[dataset]['ALL'])
        # print len(datas[dataset]['ALL']), tot
        report_ld_out.writelines('\t'.join([dataset, str(format(datas[dataset]['5']/1000000., '0,.1f'))+'M ('+str(datas[dataset]['5'] * 100/tot)+'%)', str(format(datas[dataset]['1']/1000000., '0,.1f'))+'M ('+str(datas[dataset]['1'] * 100/tot)+'%)', str(format(datas[dataset]['1_5']/1000000., '0,.1f'))+'M ('+str(datas[dataset]['1_5'] * 100/tot)+'%)', str(format(tot, '0,.0f'))])+'\n')
    report_ld_out.close()



def filter_by_maf(maf_file, region='', maf_cutoff=0.05):
    '''
    Receive a file with maf last column, return list of SNP that pass given maf
    col1 = chromosome
    col2 = position
    col3 = frequency
    '''
    POS = []
    if region != '':
        region = region.split('-')
        from_ = int(float(region[0]))+1
        to_ = int(float(region[1]))+1
    for line in open(maf_file):
        if "CHROM" not in line and "POS" not in line:
            line = line.split()
            maf = float(line[-1])
            pos = int(line[1])
            chrm = int(line[0])
            id = str(chrm)+":"+str(pos)
            if maf >= maf_cutoff:
                if region != '':
                    if from_ < pos and pos < to_:
                        POS.append(id)
                else:
                    POS.append(id)
    return POS


def report_coverage_chunk(dataset, tagNames, chunks, infoCutoff, outfile, info_format, pop_maf_format, ref_maf_format):
    """
    files:

    """
    print "Writing report per chunk to "+outfile+"..."
    output = open(outfile, 'w')
    output.writelines('\t'.join(['ID', 'CHRM', 'CHIP', 'CHUNK', 'TOTAL_IMPUTED', 'WELL_IMPUTED_NOT_CHIP', 'COMMON_IN_REF', 'COMMON_POP', 'COVERAGE', 'MIN_r2','AVERAGE_r2', 'AVERAGE_CONC'])+'\n')
    datas = {}
    dataset_data = {}
    dataset_data_ = []
    tagNames = tagNames.strip().split('==')
    chunks = chunks.strip().split('==')
    for tagName in tagNames:
        for chunk in chunks:
            # data = dat.split("==")
            # tagName = data[0]
            chunk1 = chunk.split('_')
            chrm = chunk1[0].strip()
            chunk_start = chunk1[1].strip()
            chunk_end = chunk1[2].strip()
            chunk_dataset = tagName+'__'+chrm+'__'+chunk_start+'__'+chunk_end
            info_file = info_format.format(tagName=tagName, chunk_start=chunk_start, chunk_end=chunk_end)
            pop_maf = pop_maf_format.format(tagName=tagName, chunk_start=chunk_start, chunk_end=chunk_end)
            ref_maf = ref_maf_format.format(tagName=tagName, chunk_start=chunk_start, chunk_end=chunk_end)
            ref_common_snp = filter_by_maf(ref_maf, region='', maf_cutoff=0.05) # common SNPs in the Ref panel chunk
            pop_common_snp = filter_by_maf(pop_maf, region='', maf_cutoff=0.05) # common SNPs on the chip
            info_file_data = [info.strip().split(' ') for info in open(info_file).readlines()[1:]]
            info_file_data_1 = [float(info[6]) for info in info_file_data if info[-2] == '-1']
            concord_file_data = [float(conc[-2]) for conc in info_file_data if conc[-2] != '-1']
            well_imputed = filter_info(dataset+"=="+info_file, infoCutoff, outWell_imputed=info_file+"_well_imputed", outSNP_acc=info_file+"accuracy")[1]
            well_imputed_not_chip = [chrm+":"+pos[2] for pos in well_imputed if chrm+":"+pos[2] not in pop_common_snp]
            # well_imputed_not_chip = [chrm+":"+pos[2] for pos in well_imputed]
            # imputed_not_on_chip = [pos for pos in well_imputed]
            if len(ref_common_snp) == 0 or len(pop_common_snp) == 0:
                GC = 0
            else:
                ## Coverage
                # G = 2631
                G = 6809
                L = len(well_imputed_not_chip)
                R = len(ref_common_snp)
                T = len(pop_common_snp)
                GC = round(((L/(R-T))*(G-T)+T)/float(G), 2)
                # GC1 = round((L/(R-T)) + (1-(L/(R-T))) * (T/G),2)
            ## Average r2
            if len(info_file_data_1) == 0:
                imputed_info_average = str(0)
                imputed_info_min = str(0)
            else:
                imputed_info_average = str(round(sum(info_file_data_1)/len(info_file_data_1), 2))
                imputed_info_min = str(min(info_file_data_1))
            # Concordance
            if len(concord_file_data) == 0:
                concord_file_data_average = str(0)
            else:
                concord_file_data_average = str(round(sum(concord_file_data)/len(concord_file_data), 2))
            # datas[chunk] = [str(len(info_file_data)), str(len(well_imputed))]
            dataset_data[chunk_dataset] = [tagName, chunk, str(len(well_imputed)), str(len(well_imputed_not_chip)), str(len(ref_common_snp)), str(len(pop_common_snp)), str(GC), imputed_info_min, imputed_info_average, concord_file_data_average]
            dataset_data_.append(chunk_dataset.split('__'))
    dataset_data_ = sorted(dataset_data_, key=(lambda x: (x[0], x[1], x[2])))
    dataset_data_1 = sorted([x.split('__') for x in set([ '__'.join(x[1:]) for x in dataset_data_ ])], key=(lambda x: (x[0], x[1], x[2])))
    # print dataset_data_1
    chrm_ = {}
    dataset_data_2 = {}
    for dat in dataset_data_1:
        chrm = dat[0]
        if chrm not in chrm_:
            chrm_[chrm] = 1
        else:
            chrm_[chrm] += 1
        dataset_data_2['__'.join(dat)] = '__'.join([chrm, str(chrm_[chrm])])
    print dataset_data_, dataset_data_1, dataset_data_2
    for dataset in dataset_data_:
        id = '__'.join(dataset)
        id1 = '__'.join(dataset[1:])
        output.writelines('\t'.join([dataset_data_2[id1], dataset_data_2[id1].split('_')[0]] + dataset_data[id])+'\n')
    output.close()


def report_performance_all_by_chrm(dataset, infos, infoCutoff, outfile):
    """
    files:

    """
    print "Writing report per chunk to "+outfile+"..."
    output = open(outfile, 'w')
    output.writelines('\t'.join(['CHIP', 'CHROM', 'TOTAL_IMPUTED', 'WELL_IMPUTED_NOT_CHIP', 'MIN_r2','AVERAGE_r2', 'AVERAGE_CONC'])+'\n')
    datas = {}
    for dat in infos.strip().split(','):
        data = dat.split("==")
        tagName = data[0]
        chrm = data[1]
        info_file = data[2]
        info_file_data = [info.strip().split(' ') for info in open(info_file).readlines()[1:]]
        info_file_data_1 = [float(info[6]) for info in info_file_data if info[-2] == '-1']
        concord_file_data = [float(conc[-2]) for conc in info_file_data if conc[-2] != '-1']
        well_imputed = filter_info(dataset+"=="+info_file, infoCutoff, outWell_imputed=info_file+"_well_imputed", outSNP_acc=info_file+"accuracy")[1]
        well_imputed_not_chip = [ chrm+":"+pos[2] for pos in well_imputed ]
        ## Average r2
        if len(info_file_data_1) == 0:
            imputed_info_average = str(0)
            imputed_info_min = str(0)
        else:
            imputed_info_average = str(round(sum(info_file_data_1)/len(info_file_data_1), 2))
            imputed_info_min = str(min(info_file_data_1))
        # Concordance
        if len(concord_file_data) == 0:
            concord_file_data_average = str(0)
        else:
            concord_file_data_average = str(round(sum(concord_file_data)/len(concord_file_data), 2))
        # datas[chunk] = [str(len(info_file_data)), str(len(well_imputed))]
        output.writelines('\t'.join([tagName, chrm, str(len(info_file_data)), str(len(well_imputed_not_chip)), imputed_info_min, imputed_info_average, concord_file_data_average])+'\n')
    output.close()


if args.infoFiles and args.infoCutoff:
    filter_info(args.infoFiles, args.infoCutoff, args.outWell_imputed, args.outSNP_acc)
if args.inWell_imputed and args.outWell_imputed:
    well_imputed_by_maf(args.inWell_imputed, args.outWell_imputed)
if args.inSNP_acc and args.outSNP_acc:
    acc_by_maf(args.inSNP_acc, args.outSNP_acc)
if args.report_ld and args.ldFiles:
    ld_by_maf(args.ldFiles, args.report_ld, args.inWell_imputed, args.infoCutoff)
if args.cov_chunk:
    report_coverage_chunk(args.dataset, args.tagNames, args.chunks, args.infoCutoff, args.cov_chunk, args.info_format, args.pop_maf_format, args.ref_maf_format)
if args.cov_chrm:
    report_performance_all_by_chrm(args.dataset, args.infos, args.infoCutoff, args.cov_chrm)
if args.maf_file:
    print filter_by_maf(args.maf_file, region='', maf_cutoff=0.05)
