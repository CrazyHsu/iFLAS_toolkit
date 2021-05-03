#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: identifyDEhybrid.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-04 20:24:25
Last modified: 2019-12-04 20:24:26
'''

import argparse, datetime
from multiprocessing import Pool
from itertools import combinations

from Config import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

##############################################
# def deAnalysis(comp, ngsCondDict, expFile, ngsCondFile, sampleList, refStrain):
#     if len(ngsCondDict[comp[0]]) >= 2 and len(ngsCondDict[comp[1]]) >= 2:
#         cmd = "DESeq2.R {} {} {} {} {} {}".format(expFile, ngsCondFile, comp[0], comp[1], refStrain, " ".join(sampleList))
#     else:
#         cmd = "EBseq.R {} {} {} {} {}".format(expFile, ngsCondFile, comp[0], comp[1], refStrain)
#     subprocess.call(cmd, shell=True, executable="/bin/bash")

def deAnalysis2(expFile, newColNames, tool="deseq2", selectColNames=None, conditions=None):
    if tool == "deseq2":
        cmd = "DESeq2.R {} {} {} {}".format(expFile, newColNames, selectColNames, conditions)
    else:
        cmd = "EBseq.R {} {} {} {}".format(expFile, newColNames, selectColNames, conditions)
    subprocess.call(cmd, shell=True, executable="/bin/bash")

# def identifyDEhybrid1(expFileDict=None, ngsCondFile=None, projectSample=None, analysisType="deg", dirSpec=None):
#     print str(datetime.datetime.now()) + " Identify differential expressed {}...".format("genes" if analysisType == "deg" else "transcripts")
#     prevDir = os.getcwd()
#     for projectName in expFileDict:
#         expFile = expFileDict[projectName]
#         sampleList = projectSample[projectName]
#         if expFile == None:
#             print "No expression file have been found, maybe you have provided PacBio reads without short reads, and then PacBio reads only can't be quantified! No differential expressed genes can be identified!"
#         else:
#             ngsCondDict = getSampleCond(ngsCondFile)
#             if len(ngsCondDict[projectName]) == 1:
#                 print "Only one condition have been provided, so no differential analysis will be carried out!"
#             else:
#                 deDir = os.path.join(dirSpec.out_dir, projectName, "hybrid", analysisType)
#                 resolveDir(deDir)
#                 compGroups = ngsCondDict[projectName]
#                 try:
#                     pool = Pool(processes=len(list(combinations(compGroups, 2))))
#                     multiResults = []
#                     for comp in list(combinations(compGroups, 2)):
#                         singleRes = pool.apply_async(deAnalysis, (comp, ngsCondDict, projectName, expFile, ngsCondFile, sampleList))
#                         multiResults.append(singleRes)
#                     for j in multiResults:
#                         j.wait()
#                 except Exception as e:
#                     print e
#     os.chdir(prevDir)
#     print str(datetime.datetime.now()) + " Differential expression analysis done!"

def identifyDEhybrid2(expFile, sample2bamDict, samples, projectName, dirSpec, analysisType=None):
    print str(datetime.datetime.now()) + " Identify differential expressed {} for project {}...".format(analysisType, projectName)
    prevDir = os.getcwd()

    deDir = os.path.join(dirSpec.out_dir, projectName, "hybrid", "deAnalysis", analysisType)
    resolveDir(deDir)
    newColName = [sample2bamDict["bam2sample"][i] for i in sample2bamDict["bamList"]]
    for comp1, comp2 in combinations(samples, 2):
        compList = []
        condList = []
        comp1repeats = len(comp1.ngsLeftReads.split(";"))
        comp2repeats = len(comp2.ngsLeftReads.split(";"))
        for i in range(comp1repeats):
            repeatName = "repeat" + str(i)
            compList.append(comp1.sampleName + "_" + repeatName)
            condList.append(comp1.sampleName)
        for i in range(comp2repeats):
            repeatName = "repeat" + str(i)
            compList.append(comp2.sampleName + "_" + repeatName)
            condList.append(comp2.sampleName)
        compStr = ",".join(compList)
        condStr = ",".join(condList)
        if comp1repeats == 1 and comp2repeats == 1:
            deAnalysis2(expFile, ",".join(newColName), tool="ebseq2", selectColNames=compStr, conditions=condStr)
        else:
            deAnalysis2(expFile, ",".join(newColName), tool="deseq2", selectColNames=compStr, conditions=condStr)

    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Differential expression analysis done!"

# def identifyDEhybrid(expFile, ngsCondFile, refStrain, projectName, dirSpec, compCondsDict=None, sampleNameList=None, analysisType=None):
#     print str(datetime.datetime.now()) + " Identify differential expressed {} for project {}...".format(analysisType, projectName)
#     prevDir = os.getcwd()
#
#     deDir = os.path.join(dirSpec.out_dir, projectName, "hybrid", "deAnalysis", analysisType)
#     resolveDir(deDir)
#     conditions = compCondsDict.keys()
#     try:
#         pool = Pool(processes=len(list(combinations(conditions, 2))))
#         for comp in list(combinations(conditions, 2)):
#             pool.apply_async(deAnalysis, (comp, compCondsDict, expFile, ngsCondFile, sampleNameList, refStrain))
#         pool.close()
#         pool.join()
#     except Exception as e:
#         print e
#     os.chdir(prevDir)
#     print str(datetime.datetime.now()) + " Differential expression analysis done!"

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    ngsCondFile = researchCfg.ngs_sample_cond
    projectSampleDict = {}
    expGeneFileDict = {}
    expTransFileDict = {}
    for projectName in projectsDict:
        sampleNameList = []
        for project in projects:
            if project.projectName != projectName: continue
            sampleNameList.append(project.sampleName)
        expGeneFileDict[projectName] = os.path.join(refParams.out_dir, projectName, "hybrid", "quant", "featureCounts.gene.txt")
        expTransFileDict[projectName] = os.path.join(refParams.out_dir, projectName, "hybrid", "quant", "featureCounts.transcript.txt")
        projectSampleDict[projectName] = sampleNameList
    identifyDEhybrid(expFileDict=expGeneFileDict, ngsCondFile=ngsCondFile, projectSample=projectSampleDict, refParams=refParams, analysisType="deg")
    identifyDEhybrid(expFileDict=expTransFileDict, ngsCondFile=ngsCondFile, projectSample=projectSampleDict, refParams=refParams, analysisType="det")

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)
