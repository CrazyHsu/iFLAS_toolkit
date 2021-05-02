#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: iflas.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 14:19:36
Last modified: 2021-04-29 14:19:36
'''

import sys, argparse, time
# from commonFuncs import *
from Config import *

class Logger(object):
    def __init__(self, fn):
        self.log_file = open(fn, 'w')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        if isinstance(msg,list):
            msg = ' '.join(msg)
        else:
            msg = ('{}: '+msg).format(time.ctime())
        print >> self.log_file, msg
        print msg

def iFLAS(args):
    defaultCfg = Config(args.default_cfg)
    dataToProcess = defaultCfg.dataToProcess
    refInfoParams = defaultCfg.refInfoParams
    ccsParams = defaultCfg.ccsParams
    minimap2Params = defaultCfg.minimap2Params
    optionTools = defaultCfg.optionTools
    dirSpec = defaultCfg.dirParams
    for refStrain in refInfoParams:
        refParams = refInfoParams[refStrain]
        initRefSetting(refParams=refParams, dirSpec=dirSpec)
    initSysResourceSetting(optionTools=optionTools)
    for dataObj in dataToProcess:
        refParams = refInfoParams[dataObj.ref_strain]
        dataObj.single_run_threads = len(optionTools.threads/float(len(dataToProcess)))
        if args.command == 'preproc':
            log = Logger(args.o + '.preproc.log.txt')
            log.log(sys.argv)
            from preprocess import preprocess
            preprocess(dataObj=dataObj, ccsParams=ccsParams, dirSpec=dirSpec, threads=dataObj.single_run_threads)
        if args.command == 'mapping':
            log = Logger(args.o + '.mapping.log.txt')
            log.log(sys.argv)
            from mapping import mapping
            mapping(dataObj=dataObj, minimap2Params=minimap2Params, refParams=refParams, dirSpec=dirSpec, threads=dataObj.single_run_threads)
        if args.command == 'filter':
            log = Logger(args.o + '.filter.log.txt')
            log.log(sys.argv)
            from filter import filter
            filter(args, log)
        if args.command == 'collapse':
            log = Logger(args.o + '.collapse.log.txt')
            log.log(sys.argv)
            from collapse import collapse
            collapse(args, log)
        if args.command == 'find_as':
            log = Logger(args.o + '.find_as.log.txt')
            log.log(sys.argv)
            from find_as import find_as
            find_as(args, log)
        if args.command == 'visual_as':
            log = Logger('.'.join(args.phe.split('/')[-1].split('.')[:-1]) + '.visual_as.log.txt')
            log.log(sys.argv)
            from visual_as import visual_as
            visual_as(args, log)
        if args.command == 'rank_as':
            log = Logger(args.o + '.find_as.log.txt')
            log.log(sys.argv)
            from rank_as import rank_as
            rank_as(args, log)
        if args.command == 'allelic_as':
            log = Logger(args.o + '.find_as.log.txt')
            log.log(sys.argv)
            from allelic_as import allelic_as
            allelic_as(args, log)
        if args.command == 'palen_as':
            log = Logger(args.o + '.find_as.log.txt')
            log.log(sys.argv)
            from palen_as import palen_as
            palen_as(args, log)
        if args.command == 'diff_as':
            log = Logger(args.o + '.find_as.log.txt')
            log.log(sys.argv)
            from diff_as import diff_as
            diff_as(args, log)
        if args.command == 'go':
            log = Logger(args.o + '.find_as.log.txt')
            log.log(sys.argv)
            from go import go
            go(args, log)
        if args.command == 'report':
            log = Logger(args.o + '.find_as.log.txt')
            log.log(sys.argv)
            from report import report
            report(args, log)

if __name__ == "__main__":
    USAGE = ' iFLAS: integrated Full Length Alternative Splicing analysis '

    parser = argparse.ArgumentParser(usage='%(prog)s command [options]', description=USAGE)
    subparsers = parser.add_subparsers(title='command', metavar='', dest='command', prog=parser.prog)
    parser_preprocess = subparsers.add_parser('preproc', help='Pre-process the raw PacBio/NanoPore/NGS data. When TGS and NGS data both are provide, This step will use fmlrc2 to correct the TGS read with the information in NGS', usage='%(prog)s [options]')
    # parser_preprocess.add_argument('-g', metavar='', help='Genotype file used for genome-wide genotype dimensionality reduction analysis, support plink bed format and hapmap format，the plink bed format file only needs to provide the file prefix, the hapmap format needs to be used together with the parameter -convert')
    # parser_preprocess.add_argument('-w', default=1000000, type=float, metavar='[default:1000000]', help='Window size for genome-wide genotype dimensionality reduction analysis(bp), default 1000000')
    # parser_preprocess.add_argument('-s', default=500000, type=float, metavar='[default:500000]', help='Step size for genome-wide genotype dimensionality reduction analysis(bp), default 500000')
    # parser_preprocess.add_argument('-genome_cluster', action='store_true', help='Perform genome-wide genotype dimensionality reduction analysis on genotype file')
    # parser_preprocess.add_argument('-convert', action='store_true', help='Convert hapmap genotype format to plink bed format genotype')
    # parser_preprocess.add_argument('-clump', action='store_true', help='Clumping analysis for genotype file, used to keep a subset of SNPs that are nearly uncorrelated with each other, thereby reducing the number of total SNPs to speed up genome-wide genotype dimensionality reduction analysis')
    # parser_preprocess.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads for genome cluster analysis,max threads is total number of genome chromosomes')
    # parser_preprocess.add_argument('-o', default='MODAS_genoidx_output', metavar='[default:MODAS_genoidx_output]', help='The prefix of the output file, default MODAS_genoidx_output')

    parser_mapping = subparsers.add_parser('mapping', help='Mapping the TGS/NGS reads to the reference genome with minimap2', usage='%(prog)s [options]')
    # parser_mapping.add_argument('-phe', metavar='', help='Phenotype file for phenotype preprocessing and transformation, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively')
    # parser_mapping.add_argument('-g', metavar='', help='Genotype file in plink bed format for principal component analysis, used with the parameter -pca to correct differences in phenotype caused by population structure')
    # parser_mapping.add_argument('-r', default=None, type=float, metavar='', help='Phenotype missing ratio filter cutoff')
    # parser_mapping.add_argument('-v', default=None, type=float, metavar='', help='Phenotype value filter cutoff')
    # parser_mapping.add_argument('-log2', action='store_true', help=' log2 transformation of phenotype')
    # parser_mapping.add_argument('-log10', action='store_true', help='log10 trasform phenotype')
    # parser_mapping.add_argument('-ln', action='store_true', help='ln transform phenotype')
    # parser_mapping.add_argument('-norm', action='store_true', help='Boxcox transformation of phenotype')
    # parser_mapping.add_argument('-qqnorm', action='store_true', help='Normal quantile transformation of phenotype')
    # parser_mapping.add_argument('-pca', action='store_true', help='Correct the differences in phenotype caused by population structure through PCA')
    # parser_mapping.add_argument('-o', default='MODAS_phenorm_output', metavar='[default:MODAS_phenorm_output]', help='output file prefix')

    parser_filter = subparsers.add_parser('filter', help='Screen out the low quality TGS read with junction information provided by NGS data', usage='%(prog)s [options]')
    # parser_prescreen.add_argument('-g', metavar='', help='Genotype file in plink bed format, used to calculate the kinship matrix for GWAS analysis')
    # parser_prescreen.add_argument('-genome_cluster', metavar='', help='Pseudo-genotype file for phenotype pre-screening, generated by subcommand genoidx')
    # parser_prescreen.add_argument('-phe', metavar='', help='Phenotype file for phenotype pre-screening analysis, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively')
    # parser_prescreen.add_argument('-lm_suggest_pvalue', default=None, type=float, metavar='[default:1/ genome cluster file converted pseudo-SNP number]', help='Linear model suggest pvalue for pseudo-SNP filter used for association analysis in mixed linear model, default 1/ pseudo-genotype file pseudo-snp number' )
    # parser_prescreen.add_argument('-lmm_suggest_pvalue', default=1e-5, type=float, metavar='[default:1e-5]', help=' LMM suggest pvalue for phenotype prescreen, used to obtain candidate phenotypes with significant association analysis signals and significant association regions of phenotypes, default 1e-5')
    # parser_prescreen.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads for analysis in prescreen sub command')
    # parser_prescreen.add_argument('-o', default='MODAS_prescreen_output', metavar='[default:MODAS_prescreen_output]', help='The prefix of the output file, default MODAS_prescreen_output')
    
    parser_collapse = subparsers.add_parser('collapse', help='Collapse corrected reads into high-confidence isoforms', usage='%(prog)s [options]')
    # parser_collapse.add_argument('-g', metavar='', help='Genotype file in plink bed format, used to calculate the kinship matrix for regional gwas analysis and extract snp in significant association regions of phenotypes')
    # parser_collapse.add_argument('-phe', metavar='', help='Candidate phenotype file generated by subcommand prescreen')
    # parser_collapse.add_argument('-phe_sig_qtl', metavar='', help='Significant association regions of candidate phenotypes file generated by subcommand prescreen')
    # parser_collapse.add_argument('-p1', default=1e-7, type=float, metavar='[default:1e-7]', help='Significance threshold for index SNPs, used to determine the phenotype with QTL and index snps of phenotype')
    # parser_collapse.add_argument('-p2', default=1e-6, type=float, metavar='[default:1e-6]', help='Secondary significance threshold for clumped SNPs, Used to obtain the secondary significant SNP linked to index snp to determine the reliability of the QTL. MODAS outputs the QTL containing 10 secondary significant SNPs as the phenotypic QTL result')
    # parser_collapse.add_argument('-cluster', action='store_true', help='Cluster phenotype by QTL positon')
    # parser_collapse.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads for regional gwas analysis, default 1')
    # parser_collapse.add_argument('-o', default='MODAS_regiongwas_out', metavar='[default:MODAS_regiongwas_out]', help='The prefix of the output file, default MODAS_regiongwas_out')

    parser_findAS = subparsers.add_parser('find_as', help='Identify alternative splicing(AS) type from high-confidence isoforms. Four common AS type are included: intron retention, exon skipping, alternative 3 end splicing and alternative 5 end splicing', usage='%(prog)s [options]')
    # parser_findAS.add_argument('-g', metavar='', help='Genotype file in plink bed format')
    # parser_findAS.add_argument('-exposure', metavar='', help='Exposure phenotype file, such as mTrait phenotype, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively')
    # parser_findAS.add_argument('-outcome', metavar='', help='Outcome phenotype file, such as pTrait phenotype, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively')
    # parser_findAS.add_argument('-qtl', metavar='', help='Exposure phenotype QTL file, generated by subcommand regiongwas')
    # parser_findAS.add_argument('-lm', action='store_true', help='Perform Mendelian Randomization through linear model')
    # parser_findAS.add_argument('-mlm', action='store_true', help='Perform Mendelian Randomization through mixed linear model')
    # parser_findAS.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads for Mendelian Randomization analysis')
    # parser_findAS.add_argument('-pvalue', default=1, type=float, metavar='[default:1]', help='The pvalue cutoff of Mendelian randomization analysis result output,default 1')
    # parser_findAS.add_argument('-net', action='store_true', help='Mendelian Randomization network analysis, used to identify functional modules')
    # parser_findAS.add_argument('-o', default='MODAS_mr_out', metavar='[default:MODAS_mr_out]', help='The prefix of the output file, default MODAS_mr_out')

    parser_visualAS = subparsers.add_parser('visual_as', help='Visualize the specific gene structure with details including isoform mapping, short reads coverage and AS types identified', usage='%(prog)s [options]')
    # parser_visual.add_argument('-g', metavar='', help='Genotype file in plink bed format for whole genome-wide association analysis')
    # parser_visual.add_argument('-phe', metavar='', help='Phenotype file for GWAS analysis and visualization')
    # parser_visual.add_argument('-qtl', metavar='', help='Phenotype QTL file generated by subcommand regiongwas')
    # parser_visual.add_argument('-visual', action='store_true', help='Generate web-based GWAS visualization results')
    # parser_visual.add_argument('-anno', metavar='', help='Gene annotation file for QTL annotation, used to display QTL information in visualized results')
    # parser_visual.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads genome-wide association analysis and visualization')
    # parser_visual.add_argument('-o', default='MODAS_visual_out',metavar='[default:MODAS_visual_out]', help='The prefix of the output file,default MODAS_visual_out')

    parser_rankAS = subparsers.add_parser('rank_as', help='Score the isoform by the produce of each inclusion/exclusion ratio in that isoform, and rank all the isoforms from high to low', usage='%(prog)s [options]')

    parser_allelicAS = subparsers.add_parser('allelic_as', help='Identify allelic-related AS', usage='%(prog)s [options]')

    parser_palenAS = subparsers.add_parser('palen_as', help='Identify functional poly(A) tail length related to AS', usage='%(prog)s [options]')

    parser_diffAS = subparsers.add_parser('diff_as', help='Carry out differential AS ananlysis among conditions', usage='%(prog)s [options]')

    parser_goAS = subparsers.add_parser('go', help='Perform GO enrichment analysis and plot results for the specified gene set or multiple gene sets', usage='%(prog)s [options]')

    parser_report = subparsers.add_parser('report', help='Automatic detect the plots generated in each step, and merge them into a report file', usage='%(prog)s [options]')

    if len(sys.argv) <= 2:
        sys.argv.append('-h')
    args = parser.parse_args()
    iFLAS(args)
