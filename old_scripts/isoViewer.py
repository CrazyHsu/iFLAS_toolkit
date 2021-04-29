#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: isoViewer.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-11-10 15:58:40
Last modified: 2020-01-06 20:21:03
'''

import os, re, datetime, psutil, shutil
import subprocess, argparse
import numpy as np
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool
from Bio import SeqIO

from commonFuncs import *
from commonObjs import *
from findAS import *
from Config import *
from isoViewerFuncs import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

##############################################
def isoViewer(refParams=None, tgsSample=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Visualize the gene structure compared to the reference genome for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    gpeObj = GenePredObj(refParams.ref_gpe, False)
    resolveDir(os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "isoViewer"))
    tgsIsoFile = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "ASE", "isoformGrouped.AS.confident.bed12+")
    readsDict = {}
    with open(tgsIsoFile) as f:
        for line in f:
            readStruc = ReadLineStruc(line)
            if readStruc.chrom not in readsDict:
                gene2reads = Gene2Reads(readStruc.geneName)
                gene2reads.update(readStruc)
                readsDict[readStruc.chrom] = {readStruc.strand: {readStruc.geneName: gene2reads}}
            elif readStruc.strand not in readsDict[readStruc.chrom]:
                gene2reads = Gene2Reads(readStruc.geneName)
                gene2reads.update(readStruc)
                readsDict[readStruc.chrom][readStruc.strand] = {readStruc.geneName: gene2reads}
            elif readStruc.geneName not in readsDict[readStruc.chrom][readStruc.strand]:
                gene2reads = Gene2Reads(readStruc.geneName)
                gene2reads.update(readStruc)
                readsDict[readStruc.chrom][readStruc.strand].update({readStruc.geneName: gene2reads})
            else:
                readsDict[readStruc.chrom][readStruc.strand][readStruc.geneName].update(readStruc)

    poolNum = tgsSample.threads
    # multiPlotter = []
    try:
        # pool = Pool(processes=poolNum)
        # for chrom in readsDict:
        #     for strand in readsDict[chrom]:
        #         for geneName in readsDict[chrom][strand]:
        #             gene2readsObj = readsDict[chrom][strand][geneName]
        #             sampleTargetGenePickle = pickle.dumps(gene2readsObj)
        #             if gene2readsObj.geneName in gpeObj.geneName2gpeObj:
        #                 gpeTargetGenePickle = pickle.dumps(gpeObj.geneName2gpeObj[gene2readsObj.geneName])
        #                 singlePlotter = pool.apply_async(parallelPlotterAnno, (
        #                     gene2readsObj.geneName, gpeTargetGenePickle, sampleTargetGenePickle, tgsSample, dirSpec))
        #             else:
        #                 singlePlotter = pool.apply_async(parallelPlotterNovel, (
        #                     gene2readsObj.geneName, sampleTargetGenePickle, tgsSample, dirSpec))
        #             multiPlotter.append(singlePlotter)
        # for j in multiPlotter:
        #     j.wait()
        pool = Pool(processes=poolNum)
        for chrom in readsDict:
            for strand in readsDict[chrom]:
                for geneName in readsDict[chrom][strand]:
                    gene2readsObj = readsDict[chrom][strand][geneName]
                    sampleTargetGenePickle = pickle.dumps(gene2readsObj)
                    if gene2readsObj.geneName in gpeObj.geneName2gpeObj:
                        gpeTargetGenePickle = pickle.dumps(gpeObj.geneName2gpeObj[gene2readsObj.geneName])
                        pool.apply_async(parallelPlotterAnno, (gene2readsObj.geneName, gpeTargetGenePickle,
                                                               sampleTargetGenePickle, tgsSample, dirSpec))
                    else:
                        pool.apply_async(parallelPlotterNovel, (gene2readsObj.geneName, sampleTargetGenePickle,
                                                                tgsSample, dirSpec))
        pool.close()
        pool.join()
    except Exception as e:
        print e
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Visualize the gene structure compared to the reference genome for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projectsDict = projectCfg.projectsDict
    projects = projectCfg.projects
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    tgsAnalysisList = []
    poolNum = sum([len(projectsDict[i].keys()) for i in projectsDict])
    try:
        pool = MyPool(processes=poolNum)
        multiResults = []
        for project in projects:
            if project.tgsDataDir not in tgsAnalysisList:
                tgsAnalysisList.append(project.tgsDataDir)
                singleRunRes = pool.apply_async(isoViewer, (project, refParams, projectsDict))
                multiResults.append(singleRunRes)
        for j in multiResults:
            j.wait()
    except Exception as e:
        print e

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)


