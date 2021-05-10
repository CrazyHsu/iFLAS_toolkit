#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: iflas.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 14:19:36
Last modified: 2021-04-29 14:19:36
'''

import sys, argparse, time
from commonFuncs import *
from Config import *

def mergeSample(strain2data):
    sampleMergedToProcess = {}
    from commonObjs import MergedTgsSample
    for proj in strain2data:
        for ref_strain in strain2data[proj]:
            for strain in strain2data[proj][ref_strain]:
                dataObjs = strain2data[proj][ref_strain][strain]
                sampleMerged = MergedTgsSample()
                sampleMerged.project_name = proj
                sampleMerged.sample_name = "{}_all".format(strain)
                sampleMerged.ref_strain = ref_strain
                sampleMerged.strain = strain
                sampleMerged.tgs_plat = dataObjs[0].tgs_plat
                sampleMerged.strategy = dataObjs[0].strategy
                sampleMerged.data_processed_location = [x.data_processed_location for x in dataObjs]

                minLeftRepeats = min([len(x.ngs_left_reads.split(";")) for x in dataObjs])
                leftReads = []
                leftReadsRepeats = []
                for x1 in range(minLeftRepeats):
                    for x in dataObjs:
                        leftReads.append(x.ngs_left_reads.split(";")[x1])
                    leftReadsRepeats.append(",".join(leftReads))
                sampleMerged.ngs_left_reads = ";".join(leftReadsRepeats)
                rightReads = []
                rightReadsRepeats = []
                minRightRepeats = min([len(x.ngs_right_reads.split(";")) for x in dataObjs])
                for x1 in range(minRightRepeats):
                    for x in dataObjs:
                        rightReads.append(x.ngs_right_reads.split(";")[x1])
                    rightReadsRepeats.append(",".join(rightReads))
                sampleMerged.ngs_right_reads = ";".join(rightReadsRepeats)
                if not sampleMerged.ngs_left_reads and sampleMerged.ngs_right_reads or sampleMerged.ngs_left_reads and not sampleMerged.ngs_right_reads:
                    sampleMerged.ngs_reads_paired = "single"
                else:
                    sampleMerged.ngs_reads_paired = "paired"

                sampleMerged.use_fmlrc2 = dataObjs[0].use_fmlrc2

                if proj not in sampleMergedToProcess:
                    sampleMergedToProcess[proj] = {ref_strain: {strain: [sampleMerged]}}
                elif ref_strain not in sampleMergedToProcess[proj]:
                    sampleMergedToProcess[proj][ref_strain] = {strain: [sampleMerged]}
                elif strain not in sampleMergedToProcess[proj][ref_strain]:
                    sampleMergedToProcess[proj][ref_strain][strain] = [sampleMerged]
                else:
                    sampleMergedToProcess[proj][ref_strain][strain].append(sampleMerged)
    return sampleMergedToProcess

def iflas(args):
    defaultCfg = Config(args.default_cfg)
    dataToProcess = defaultCfg.dataToProcess
    refInfoParams = defaultCfg.refInfoParams
    ccsParams = defaultCfg.ccsParams
    minimap2Params = defaultCfg.minimap2Params
    collapseParams = defaultCfg.collapseParams
    optionTools = defaultCfg.optionTools
    dirSpec = defaultCfg.dirParams
    # from commonFuncs import initRefSetting, initSysResourceSetting
    for refStrain in refInfoParams:
        refParams = refInfoParams[refStrain]
        initRefSetting(refParams=refParams, dirSpec=dirSpec)
    initSysResourceSetting(optionTools=optionTools)
    strain2data = {}
    for dataObj in dataToProcess:
        dataObj.single_run_threads = len(optionTools.threads/float(len(dataToProcess)))
        if args.command == 'preproc':
            from preprocess import preprocess
            preprocess(dataObj=dataObj, ccsParams=ccsParams, dirSpec=dirSpec, threads=dataObj.single_run_threads)
            if dataObj.project_name not in strain2data:
                strain2data[dataObj.project_name] = {dataObj.ref_strain: {dataObj.strain: [dataObj]}}
            elif dataObj.ref_strain not in strain2data[dataObj.project_name]:
                strain2data[dataObj.project_name][dataObj.ref_strain] = {dataObj.strain: [dataObj]}
            elif dataObj.strain not in strain2data[dataObj.project_name][dataObj.ref_strain]:
                strain2data[dataObj.project_name][dataObj.ref_strain][dataObj.strain] = [dataObj]
            else:
                strain2data[dataObj.project_name][dataObj.ref_strain][dataObj.strain].append(dataObj)
    if optionTools.merge_data_from_same_strain:
        sampleMergedToProcess = mergeSample(strain2data)
        for proj in sampleMergedToProcess:
            for ref_strain in sampleMergedToProcess[proj]:
                for strain in sampleMergedToProcess[proj][ref_strain]:
                    for dataObj in sampleMergedToProcess[proj][ref_strain][strain]:
                        refParams = refInfoParams[ref_strain]
                        dataObj.single_run_threads = len(optionTools.threads / float(len(sampleMergedToProcess)))
                        if args.command == 'mapping':
                            from mapping import mapping
                            mapping(dataObj=dataObj, minimap2Params=minimap2Params, refParams=refParams,
                                    dirSpec=dirSpec, threads=dataObj.single_run_threads)
                        if args.command == 'filter':
                            from filter import filter
                            filter(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                        if args.command == 'collapse':
                            from collapse import collapse
                            collapse(dataObj=dataObj, collapseParams=collapseParams, refParams=refParams,
                                     dirSpec=dirSpec, threads=dataObj.single_run_threads)
                        if args.command == 'find_as':
                            from find_as import find_as
                            find_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                            from find_pa import find_pa
                            find_pa(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                            from charaterize_as import charaterize_as
                            charaterize_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                        if args.command == 'visual_as':
                            from visual_as import visual_as
                            targetGenes = args.genes
                            visual_as(dataObj=dataObj, targetGenes=targetGenes, refParams=refParams, dirSpec=dirSpec)
                        if args.command == 'rank_as':
                            from rank_as import rank_as
                            rank_as(dataObj=dataObj, dirSpec=dirSpec)
                        if args.command == 'allelic_as':
                            from allelic_as import allelic_as
                            allelic_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                        if args.command == 'palen_as':
                            from palen_as import palen_as
                            palen_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec,
                                     sampleMerged=optionTools.merge_data_from_same_strain)
    else:
        for dataObj in dataToProcess:
            refParams = refInfoParams[dataObj.ref_strain]
            dataObj.single_run_threads = len(optionTools.threads / float(len(dataToProcess)))
            if args.command == 'mapping':
                from mapping import mapping
                mapping(dataObj=dataObj, minimap2Params=minimap2Params, refParams=refParams, dirSpec=dirSpec, threads=dataObj.single_run_threads)
            if args.command == 'filter':
                from filter import filter
                filter(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
            if args.command == 'collapse':
                from collapse import collapse
                collapse(dataObj=dataObj, collapseParams=collapseParams, refParams=refParams, dirSpec=dirSpec, threads=dataObj.single_run_threads)
            if args.command == 'find_as':
                # from find_charaterize_as_functions import *
                from find_as import find_as
                find_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                from find_pa import find_pa
                find_pa(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                from charaterize_as import charaterize_as
                charaterize_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
            if args.command == 'visual_as':
                from visual_as import visual_as
                targetGenes = args.genes
                visual_as(dataObj=dataObj, targetGenes=targetGenes, refParams=refParams, dirSpec=dirSpec)
            if args.command == 'rank_as':
                from rank_as import rank_as
                rank_as(dataObj=dataObj, dirSpec=dirSpec)
            if args.command == 'allelic_as':
                from allelic_as import allelic_as
                allelic_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
            if args.command == 'palen_as':
                from palen_as import palen_as
                palen_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, sampleMerged=optionTools.merge_data_from_same_strain)
            if args.command == 'report':
                from report import report
                report(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
        if args.command == 'diff_as':
            from diff_as import diff_as, diff_as1
            # diff_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
            compCond = args.compCond
            diff_as1(dataToProcess, compCondFile=compCond, dirSpec=dirSpec, sampleMerged=optionTools.merge_data_from_same_strain)
        if args.command == 'go':
            from go import go
            go(args)


if __name__ == "__main__":
    USAGE = ' iFLAS: integrated Full Length Alternative Splicing analysis '

    parser = argparse.ArgumentParser(usage='%(prog)s command [options]', description=USAGE)
    subparsers = parser.add_subparsers(title='command', metavar='', dest='command', prog=parser.prog)
    parser_preprocess = subparsers.add_parser('preproc', help='Pre-process the raw PacBio/NanoPore/NGS data. When TGS and NGS data both are provide, This step will use fmlrc2 to correct the TGS read with the information in NGS', usage='%(prog)s [options]')
    parser_preprocess.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    parser_mapping = subparsers.add_parser('mapping', help='Mapping the TGS/NGS reads to the reference genome with minimap2', usage='%(prog)s [options]')
    parser_mapping.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    parser_filter = subparsers.add_parser('filter', help='Screen out the low quality TGS read with junction information provided by NGS data', usage='%(prog)s [options]')
    parser_filter.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    parser_collapse = subparsers.add_parser('collapse', help='Collapse corrected reads into high-confidence isoforms', usage='%(prog)s [options]')
    parser_collapse.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    parser_findAS = subparsers.add_parser('find_as', help='Identify alternative splicing(AS) type from high-confidence isoforms. Four common AS type are included: intron retention, exon skipping, alternative 3 end splicing and alternative 5 end splicing', usage='%(prog)s [options]')
    parser_findAS.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    parser_visualAS = subparsers.add_parser('visual_as', help='Visualize the specific gene structure with details including isoform mapping, short reads coverage and AS types identified', usage='%(prog)s [options]')
    parser_visualAS.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")
    parser_visualAS.add_argument('-g', dest="genes", help="The gene list separated by comma or a single file contain genes one per line.")

    parser_rankAS = subparsers.add_parser('rank_as', help='Score the isoform by the produce of each inclusion/exclusion ratio in that isoform, and rank all the isoforms from high to low', usage='%(prog)s [options]')
    parser_rankAS.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    parser_allelicAS = subparsers.add_parser('allelic_as', help='Identify allelic-related AS', usage='%(prog)s [options]')
    parser_allelicAS.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    parser_palenAS = subparsers.add_parser('palen_as', help='Identify functional poly(A) tail length related to AS', usage='%(prog)s [options]')
    parser_palenAS.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    parser_diffAS = subparsers.add_parser('diff_as', help='Carry out differential AS ananlysis among conditions', usage='%(prog)s [options]')
    parser_diffAS.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")
    parser_diffAS.add_argument('-d', dest="compCond", help="The condition file used to detect differential AS between samples.")

    parser_goAS = subparsers.add_parser('go', help='Perform GO enrichment analysis and plot results for the specified gene set or multiple gene sets', usage='%(prog)s [options]')
    parser_goAS.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")
    parser_goAS.add_argument('-tg', dest="targetGeneFile", help="The target gene file or file list separated by comma.")
    parser_goAS.add_argument('-bg', dest="gene2goFile", help="The mapping file between gene and go term.")
    parser_goAS.add_argument('-s', dest="sampleName", help="The sample name used plot the track, multi-sample should be separated by commma.")

    parser_report = subparsers.add_parser('report', help='Automatic detect the plots generated in each step, and merge them into a report file', usage='%(prog)s [options]')
    parser_report.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    if len(sys.argv) <= 2:
        sys.argv.append('-h')
    args = parser.parse_args()
    iflas(args)
