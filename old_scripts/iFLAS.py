#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: iFLAS.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-08 17:12:05
Last modified: 2019-12-08 17:12:06
'''

import argparse, datetime, psutil, subprocess, pybedtools
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool

from Config import *
from commonFuncs import NoDaemonProcess, MyPool
from retrieveTGSdata import *
from preprocessAndCorrFuncs import *
from filterAndCollapseFuncs import *
from evaluationFuncs import *
from paAnalysisFuncs import *
from identifyASEFuncs import *
from isoViewerFuncs import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-p", dest="project", type=str, default="project",
                    help="The project name")
parser.add_argument("-g", dest="group", type=str, default="group",
                    help="The group which the analysis belongs to. For single running, this parameter can be missed")
parser.add_argument("-t", dest="tgs_dir", type=str, default="tgsDir",
                    help="The directory which cantain the TGS raw reads to be analysis")
parser.add_argument("-n", dest="ngs_cfg_file", type=str, default=None,
                    help="The directory which contain the NGS raw reads to be analysis")
parser.add_argument("-plat", dest="tgs_plat", type=str, default="PacBio",
                    help="The sequence strategy to carry out TGS sequencing. Default: PacBio")
parser.add_argument("-pac_plat", dest="pac_plat", type=str, default="sequel",
                    help="Specify the strategy that used to carry out PacBio sequencing, Sequel or RSII. Default: Sequel")
parser.add_argument("-pair", dest="ngs_pair", action="store_true",
                    help="If NGS raw reads provided, specify this parameter to determine the reads of NGS is paired")
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
parser.add_argument("-cc", dest="command", action="store_true",
                    help="Use this parameter to specify whether the program runs from command line or from "
                         "config file. Note: when the sequence information specified from command line the program will "
                         "ignore the information in config file and only be executed once")
args = parser.parse_args()

##############################################
def singleRun1(refParams=None, tgsSample=None, mecatParams=None, optionTools=None, lncRNAcpatParams=None, dirSpec=None):
    # initSysSetting(refParams=refParams, projectsDict=projectsDict, threads=threads)
    prevDir = os.getcwd()
    workDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName)
    resolveDir(workDir)
    if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
        preprocessNGSdata(tgsSample=tgsSample)
        pass

    if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
        renameNGSdata2fastp(tgsSample=tgsSample)

    # from processRawReadsAndCorr import processRawReadsAndCorrection, correctionAndAlignToReference
    # processRawReadsAndCorrection(tgsSample=tgsSample, refParams=refParams, ccsParams=ccsParams, mecatParams=mecatParams,
    #                              useProovread=optionTools.use_proovread, projectsDict=projectsDict)
    from processRawReadsAndCorr import correctionAndAlignToReference
    correctionAndAlignToReference(tgsSample=tgsSample, refParams=refParams, mecatParams=mecatParams,
                                  useProovread=optionTools.use_proovread, dirSpec=dirSpec)

    if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
        processReadsFromNGS(tgsSample=tgsSample, refParams=refParams, dirSpec=dirSpec)

    # if optionTools.use_cogent:
    #     # makeBadMappedTranscriptome()
    #     pass
    #
    from filterAndCollapse import filterAndCollpase
    filterAndCollpase(refParams=refParams, tgsSample=tgsSample, optionTools=optionTools, dirSpec=dirSpec)
    #
    # if optionTools.find_fusion:
    #     from fusionTransAnalysis import fusionTransAnalysis
    #     fusionTransAnalysis(tgsSample=tgsSample, dirSpec=dirSpec)
    #
    from evaluation import evaluation
    evaluation(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)
    #
    # if lncRNAcpatParams[tgsSample.refStrain].coding_seq == None or lncRNAcpatParams[tgsSample.refStrain].noncoding_seq == None:
    #     print "No coding or noncoding information provided, lncRNA identification can't be carried out!"
    # else:
    #     from identifyLncRNA import lncrnaCPAT
    #     lncrnaCPAT(refParams=refParams, tgsSample=tgsSample, lncRNAcpat=lncRNAcpatParams, dirSpec=dirSpec)
    #
    from paAnalysis import paAnalysis
    paAnalysis(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)
    # #
    from identifyASE import identifyASE, characterizeASE
    identifyASE(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)

    from isoViewer import isoViewer
    isoViewer(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)
    #
    # if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
    #     from indepComb import indepComb
    #     indepComb(tgsSample=tgsSample, dirSpec=dirSpec)
    # #
    # # # poly-A length related AS
    # from palenRelatedAS import palenRelatedAS
    # palenRelatedAS(tgsSample=tgsSample, dirSpec=dirSpec)
    # #
    # # allele specific related AS
    # from alleleSpecific import alleleSpecific, runPhaser, makeAbundanceFile, getASpairedIsoforms
    # alleleSpecific(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)
    os.chdir(prevDir)

def integratedRun1(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refInfoParams = defaultCfg.refInfoParams
    mecatParams = defaultCfg.mecatParams
    ccsParams = defaultCfg.ccsParams
    optionTools = defaultCfg.optionTools
    lncRNAcpatParams = defaultCfg.lncRNAcpatParams
    dirSpec = defaultCfg.dirParams

    projectCfg1 = ProjectCfg1(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    tgsDict = projectCfg1.tgsDict
    # tgsList = projectCfg1.tgsList
    # tgsStrain2sample = projectCfg1.tgsStrain2sample
    tgsPlat2sample = projectCfg1.tgsPlat2sample
    # tgsRefStrain2sample = projectCfg1.tgsRefStrain2sample
    tgsGroupedSample = projectCfg1.tgsGroupedSample

    ngsDict = projectCfg1.ngsDict

    for refStrain in refInfoParams:
        refParams = refInfoParams[refStrain]
        initRefSetting(refParams=refParams, dirSpec=dirSpec)
        # initSysSetting(refParams=refParams, optionTools=optionTools, dirSpec=dirSpec)
    initResourceSetting(optionTools=optionTools)

    try:
        retrieveTGSdata(tgsPlat2sample, dirSpec, refInfoParams, ccsParams, optionTools)
    except Exception as e:
        raise e

    sampleDict = {}
    myproject2refStrain = {}
    for projectName in tgsDict:
        sampleDict[projectName] = []
        myproject2refStrain[projectName] = []
        for tgsPlat in tgsPlat2sample[projectName]:
            if tgsPlat == "pacbio":
                demultiplexDir = os.path.join(dirSpec.out_dir, projectName, "retrieveData", "pacbioDemultiplexed")
                for refStrain in tgsDict[projectName]:
                    for strainName in tgsDict[projectName][refStrain]:
                        myproject2refStrain[projectName][refStrain] = []
                        for groupedSampleName in tgsGroupedSample[projectName]:
                            for i in groupedSampleName.split("+"):
                                sampleName = strainName + "_" + i
                                ngsSample = ngsDict[projectName][strainName][i]

                                mySample = MergedTgsSample()
                                mySample.projectName = projectName
                                mySample.sampleName = sampleName
                                mySample.refStrain = refStrain
                                mySample.getUniqName()
                                mySample.tgsPlat = tgsPlat
                                mySample.tgsStrategy = tgsGroupedSample[projectName][groupedSampleName].tgsStrategy
                                mySample.tgsProcessedData = os.path.join(demultiplexDir, "{}.flnc.bam".format(mySample.uniqName))
                                mySample.tgsPrimer = os.path.join(demultiplexDir, "{}.primers.fa".format(mySample.uniqName))
                                mySample.ngsLeftReads = ngsSample.leftReads
                                mySample.ngsRightReads = ngsSample.rightReads
                                mySample.ngsPaired = ngsSample.ngsReadPair
                                mySample.ngsReadsLength = ngsSample.readLength

                                sampleDict[projectName].append(mySample)
                                myproject2refStrain[projectName][refStrain].append(mySample)
            else:
                nanoporeDir = os.path.join(dirSpec.out_dir, projectName, "retrieveData", "nanoporeRetrieve")
                for refStrain in tgsDict[projectName]:
                    for strainName in tgsDict[projectName][refStrain]:
                        myproject2refStrain[projectName][refStrain] = []
                        for groupedSampleName in tgsGroupedSample[projectName]:
                            sampleName = strainName + "_" + groupedSampleName
                            ngsSample = ngsDict[projectName][strainName][groupedSampleName]

                            mySample = MergedTgsSample()
                            mySample.projectName = projectName
                            mySample.sampleName = sampleName
                            mySample.refStrain = refStrain
                            mySample.getUniqName()
                            mySample.tgsPlat = tgsPlat
                            mySample.tgsStrategy = tgsGroupedSample[projectName][groupedSampleName].tgsStrategy
                            mySample.tgsProcessedData = os.path.join(nanoporeDir, "{}", "nanoporeRawSeq.fq".format(sampleName))
                            mySample.ngsLeftReads = ngsSample.leftReads
                            mySample.ngsRightReads = ngsSample.rightReads
                            mySample.ngsPaired = ngsSample.ngsReadPair
                            mySample.ngsReadsLength = ngsSample.readLength

                            sampleDict[projectName].append(mySample)
                            myproject2refStrain[projectName][refStrain].append(mySample)


    for projectName in sampleDict:
        for tmpRun in sampleDict[projectName]:
            refParams = refInfoParams[tmpRun.refStrain]
            if tmpRun.ngsLeftReads or tmpRun.ngsRightReads:
                makeHisat2Index(refParams=refParams, dirSpec=dirSpec, optionTools=optionTools)

    for projectName in sampleDict:
        for tmpRun in sampleDict[projectName]:
            refParams = refInfoParams[tmpRun.refStrain]
            tmpRun.threads = int(optionTools.threads / float(len(sampleDict[projectName])))
            tmpRun.memory = str(int(optionTools.memory[0:-1]) / float(len(sampleDict[projectName]))) + "M"
            preprocessNGSdata(tgsSample=tmpRun)
            singleRun1(refParams=refParams, tgsSample=tmpRun, mecatParams=mecatParams, optionTools=optionTools, lncRNAcpatParams=lncRNAcpatParams, dirSpec=dirSpec)

    ngsCondFile = researchCfg.ngs_cfg_file
    ngsCondDict = getSampleCond(ngsCondFile)
    for projectName in myproject2refStrain:
        for refStrain in myproject2refStrain[projectName]:
            refParams = refInfoParams[refStrain]
            projects = myproject2refStrain[projectName][refStrain]
            sampleNameList = [i.sampleName for i in projects]
            from quantHybrid import quantifyHybrid
            geneFeatureCountFile, transFeatureCountFile = quantifyHybrid(refParams=refParams, projects=projects, dirSpec=dirSpec, projectName=projectName, refStrain=refStrain)

            from identifyDEhybrid import identifyDEhybrid, deAnalysis
            compCondsDict = ngsCondDict[projectName][refStrain]
            identifyDEhybrid(geneFeatureCountFile, ngsCondFile, refStrain, projectName, dirSpec,
                              compCondsDict=compCondsDict, sampleNameList=sampleNameList, analysisType="genes")
            identifyDEhybrid(geneFeatureCountFile, ngsCondFile, refStrain, projectName, dirSpec,
                              compCondsDict=compCondsDict, sampleNameList=sampleNameList, analysisType="transcripts")

            from identifyDAShybrid import identifyDAShybrid1, dasAnalysis
            identifyDAShybrid1(projectName, dirSpec, compCondsDict=compCondsDict, projects=projects)

            from strainSpecificAS import getStrainSpecificComparion, specificAS
            getStrainSpecificComparion(dirSpec, projects=projects)

def test_integratedRun1(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refInfoParams = defaultCfg.refInfoParams
    mecatParams = defaultCfg.mecatParams
    ccsParams = defaultCfg.ccsParams
    optionTools = defaultCfg.optionTools
    lncRNAcpatParams = defaultCfg.lncRNAcpatParams
    dirSpec = defaultCfg.dirParams

    projectCfg1 = ProjectCfg1(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    tgsDict = projectCfg1.tgsDict
    # tgsList = projectCfg1.tgsList
    # tgsStrain2sample = projectCfg1.tgsStrain2sample
    tgsPlat2sample = projectCfg1.tgsPlat2sample
    # tgsRefStrain2sample = projectCfg1.tgsRefStrain2sample
    tgsGroupedSample = projectCfg1.tgsGroupedSample

    ngsDict = projectCfg1.ngsDict

    for refStrain in refInfoParams:
        refParams = refInfoParams[refStrain]
        initRefSetting(refParams=refParams, dirSpec=dirSpec)
        # initSysSetting(refParams=refParams, optionTools=optionTools, dirSpec=dirSpec)
    initResourceSetting(optionTools=optionTools)

    # retrieveTGSdata(tgsPlat2sample, dirSpec, refInfoParams, ccsParams, optionTools)

    sampleDict = {}
    myproject2refStrain = {}
    for projectName in tgsDict:
        sampleDict[projectName] = []
        myproject2refStrain[projectName] = {}
        for tgsPlat in tgsPlat2sample[projectName]:
            if tgsPlat == "pacbio":
                demultiplexDir = os.path.join(dirSpec.out_dir, projectName, "retrieveData", "pacbioDemultiplexed")
                for refStrain in tgsDict[projectName]:
                    myproject2refStrain[projectName][refStrain] = []
                    for strainName in tgsDict[projectName][refStrain]:
                        for groupedSampleName in tgsGroupedSample[projectName]:
                            for i in groupedSampleName.split("+"):
                                if i not in ngsDict[projectName][strainName]: continue
                                sampleName = strainName + "_" + i
                                ngsSample = ngsDict[projectName][strainName][i]

                                mySample = MergedTgsSample()
                                mySample.projectName = projectName
                                mySample.sampleName = sampleName
                                mySample.refStrain = refStrain
                                mySample.getUniqName()
                                mySample.tgsPlat = tgsPlat
                                mySample.tgsStrategy = tgsGroupedSample[projectName][groupedSampleName].tgsStrategy
                                tgsUniqName = "{}_{}_{}".format(projectName, strainName, groupedSampleName)
                                mySample.tgsProcessedData = os.path.join(demultiplexDir, "{}.{}.flnc.bam".format(tgsUniqName, i))
                                if "endo" in i:
                                    cmd = "samtools cat {}/{}.{}.flnc.bam {}/test1_All_endo_B73_endo+Ki11_endo+B73_Ki11_endo+Ki11_B73_endo.{}_endo.flnc.bam > {}/{}.{}.merged.flnc.bam".format(demultiplexDir, tgsUniqName, i, demultiplexDir, strainName, demultiplexDir, tgsUniqName, i)
                                    # subprocess.call(cmd, shell=True)
                                    mySample.tgsProcessedData = os.path.join(demultiplexDir, "{}.{}.merged.flnc.bam".format(tgsUniqName, i))
                                mySample.tgsPrimer = os.path.join(demultiplexDir, "{}.{}.primers.fa".format(tgsUniqName, i))
                                mySample.ngsLeftReads = ngsSample.leftReads
                                mySample.ngsRightReads = ngsSample.rightReads
                                mySample.ngsPaired = ngsSample.ngsReadPair
                                mySample.ngsReadsLength = ngsSample.readLength

                                sampleDict[projectName].append(mySample)
                                myproject2refStrain[projectName][refStrain].append(mySample)
            else:
                nanoporeDir = os.path.join(dirSpec.out_dir, projectName, "retrieveData", "nanoporeRetrieve")
                for refStrain in tgsDict[projectName]:
                    for strainName in tgsDict[projectName][refStrain]:
                        myproject2refStrain[projectName][refStrain] = []
                        for groupedSampleName in tgsGroupedSample[projectName]:
                            sampleName = strainName + "_" + groupedSampleName
                            ngsSample = ngsDict[projectName][strainName][groupedSampleName]

                            mySample = MergedTgsSample()
                            mySample.projectName = projectName
                            mySample.sampleName = sampleName
                            mySample.refStrain = refStrain
                            mySample.getUniqName()
                            mySample.tgsPlat = tgsPlat
                            mySample.dataLocation = tgsGroupedSample[projectName][groupedSampleName].dataLocation
                            mySample.tgsStrategy = tgsGroupedSample[projectName][groupedSampleName].tgsStrategy
                            mySample.tgsProcessedData = os.path.join(nanoporeDir, sampleName, "nanoporeRawSeq.fq")
                            mySample.ngsLeftReads = ngsSample.leftReads
                            mySample.ngsRightReads = ngsSample.rightReads
                            mySample.ngsPaired = ngsSample.ngsReadPair
                            mySample.ngsReadsLength = ngsSample.readLength

                            sampleDict[projectName].append(mySample)
                            myproject2refStrain[projectName][refStrain].append(mySample)

    for projectName in sampleDict:
        for tmpRun in sampleDict[projectName]:
            refParams = refInfoParams[tmpRun.refStrain]
            if tmpRun.ngsLeftReads or tmpRun.ngsRightReads:
                makeHisat2Index(refParams=refParams, dirSpec=dirSpec, optionTools=optionTools)

    # pybedtools.set_tempdir(dirSpec.tmp_dir)
    # projectNum = sum([len(sampleDict[i]) for i in sampleDict])
    # pool = MyPool(processes=projectNum)
    # for projectName in sampleDict:
    #     for tmpRun in sampleDict[projectName]:
    #         refParams = refInfoParams[tmpRun.refStrain]
    #         tmpRun.threads = int(optionTools.threads / float(projectNum))
    #         tmpRun.memory = str(int(optionTools.memory[0:-1]) / float(projectNum)) + "M"
    #         pool.apply_async(singleRun1, (refParams, tmpRun, mecatParams, optionTools, lncRNAcpatParams, dirSpec))
    # pool.close()
    # pool.join()


    mergedSampleList = []
    for projectName in tgsDict:
        for tgsPlat in tgsPlat2sample[projectName]:
            if tgsPlat == "pacbio":
                demultiplexDir = os.path.join(dirSpec.out_dir, projectName, "retrieveData", "pacbioDemultiplexed")
                for refStrain in tgsDict[projectName]:
                    for strainName in tgsDict[projectName][refStrain]:
                        leftReadsR1, leftReadsR2, rightReadsR1, rightReadsR2 = [], [], [], []
                        tgsData = []
                        for groupedSampleName in tgsGroupedSample[projectName]:
                            for i in groupedSampleName.split("+"):
                                if i not in ngsDict[projectName][strainName]: continue
                                sampleName = strainName + "_" + i
                                ngsSample = ngsDict[projectName][strainName][i]

                                mySample = MergedTgsSample()
                                mySample.projectName = projectName
                                mySample.sampleName = sampleName
                                mySample.refStrain = refStrain
                                mySample.getUniqName()
                                mySample.tgsPlat = tgsPlat
                                mySample.tgsStrategy = tgsGroupedSample[projectName][groupedSampleName].tgsStrategy
                                tgsUniqName = "{}_{}_{}".format(projectName, strainName, groupedSampleName)
                                mySample.tgsProcessedData = os.path.join(demultiplexDir, "{}.{}.flnc.bam".format(tgsUniqName, i))
                                mySample.tgsPrimer = os.path.join(demultiplexDir, "{}.{}.primers.fa".format(tgsUniqName, i))
                                mySample.ngsLeftReads = ngsSample.leftReads
                                mySample.ngsRightReads = ngsSample.rightReads
                                mySample.ngsPaired = ngsSample.ngsReadPair
                                mySample.ngsReadsLength = ngsSample.readLength
                                leftReadsR1.append(ngsSample.leftReads.split(";")[0])
                                leftReadsR2.append(ngsSample.leftReads.split(";")[1])
                                rightReadsR1.append(ngsSample.rightReads.split(";")[0])
                                rightReadsR2.append(ngsSample.rightReads.split(";")[1])
                                tgsData.append(os.path.join(demultiplexDir, "{}.{}.flnc.bam".format(tgsUniqName, i)))

                                # sampleDict[projectName].append(mySample)
                                # myproject2refStrain[projectName][refStrain].append(mySample)
                            allSampleMerged = MergedTgsSample()
                            if strainName == "B73":
                                allSampleMerged.projectName = projectName
                                allSampleMerged.sampleName = "B73_all"
                                allSampleMerged.refStrain = refStrain
                                allSampleMerged.getUniqName()
                                allSampleMerged.tgsPlat = tgsPlat
                                allSampleMerged.tgsStrategy = "sequel"
                                allSampleMerged.ngsLeftReads = ",".join(leftReadsR1) + ";" + ",".join(leftReadsR2)
                                allSampleMerged.ngsRightReads = ",".join(rightReadsR1) + ";" + ",".join(rightReadsR2)
                                allSampleMerged.ngsPaired = "paired"
                                allSampleMerged.ngsReadsLength = 151
                                allSampleMerged.tgsProcessedData = "{}/B73_all.flnc.bam".format(demultiplexDir)
                                cmd = "samtools cat {} {}/test1_All_endo_B73_endo+Ki11_endo+B73_Ki11_endo+Ki11_B73_endo.B73_endo.flnc.bam > {}".format(" ".join(tgsData), demultiplexDir, allSampleMerged.tgsProcessedData)
                                # subprocess.call(cmd, shell=True)
                            elif strainName == "Ki11":
                                allSampleMerged.projectName = projectName
                                allSampleMerged.sampleName = "Ki11_all"
                                allSampleMerged.refStrain = refStrain
                                allSampleMerged.getUniqName()
                                allSampleMerged.tgsPlat = tgsPlat
                                allSampleMerged.tgsStrategy = "sequel"
                                allSampleMerged.ngsLeftReads = ",".join(leftReadsR1) + ";" + ",".join(leftReadsR2)
                                allSampleMerged.ngsRightReads = ",".join(rightReadsR1) + ";" + ",".join(rightReadsR2)
                                allSampleMerged.ngsPaired = "paired"
                                allSampleMerged.ngsReadsLength = 151
                                allSampleMerged.tgsProcessedData = "{}/Ki11_all.flnc.bam".format(demultiplexDir)
                                cmd = "samtools cat {} {}/test1_All_endo_B73_endo+Ki11_endo+B73_Ki11_endo+Ki11_B73_endo.Ki11_endo.flnc.bam > {}".format(
                                    " ".join(tgsData), demultiplexDir, allSampleMerged.tgsProcessedData)
                                # subprocess.call(cmd, shell=True)
                            elif strainName == "B73_Ki11":
                                allSampleMerged.projectName = projectName
                                allSampleMerged.sampleName = "B73_Ki11_all"
                                allSampleMerged.refStrain = refStrain
                                allSampleMerged.getUniqName()
                                allSampleMerged.tgsPlat = tgsPlat
                                allSampleMerged.tgsStrategy = "sequel"
                                allSampleMerged.ngsLeftReads = ",".join(leftReadsR1) + ";" + ",".join(leftReadsR2)
                                allSampleMerged.ngsRightReads = ",".join(rightReadsR1) + ";" + ",".join(rightReadsR2)
                                allSampleMerged.ngsPaired = "paired"
                                allSampleMerged.ngsReadsLength = 151
                                allSampleMerged.tgsProcessedData = "{}/B73_Ki11_all.flnc.bam".format(demultiplexDir)
                                cmd = "samtools cat {} {}/test1_All_endo_B73_endo+Ki11_endo+B73_Ki11_endo+Ki11_B73_endo.B73_Ki11_endo.flnc.bam > {}".format(
                                    " ".join(tgsData), demultiplexDir, allSampleMerged.tgsProcessedData)
                                # subprocess.call(cmd, shell=True)
                            else:
                                allSampleMerged.projectName = projectName
                                allSampleMerged.sampleName = "Ki11_B73_all"
                                allSampleMerged.refStrain = refStrain
                                allSampleMerged.getUniqName()
                                allSampleMerged.tgsPlat = tgsPlat
                                allSampleMerged.tgsStrategy = "sequel"
                                allSampleMerged.ngsLeftReads = ",".join(leftReadsR1) + ";" + ",".join(leftReadsR2)
                                allSampleMerged.ngsRightReads = ",".join(rightReadsR1) + ";" + ",".join(rightReadsR2)
                                allSampleMerged.ngsPaired = "paired"
                                allSampleMerged.ngsReadsLength = 151
                                allSampleMerged.tgsProcessedData = "{}/Ki11_B73_all.flnc.bam".format(demultiplexDir)
                                cmd = "samtools cat {} {}/test1_All_endo_B73_endo+Ki11_endo+B73_Ki11_endo+Ki11_B73_endo.Ki11_B73_endo.flnc.bam > {}".format(
                                    " ".join(tgsData), demultiplexDir, allSampleMerged.tgsProcessedData)
                                # subprocess.call(cmd, shell=True)
                            mergedSampleList.append(allSampleMerged)

    # pool = MyPool(processes=len(mergedSampleList))
    # for i in mergedSampleList:
    #     tmpSample = i
    #     refParams = refInfoParams[tmpSample.refStrain]
    #     tmpSample.threads = int(optionTools.threads / float(len(mergedSampleList)))
    #     tmpSample.memory = str(int(optionTools.memory[0:-1]) / float(len(mergedSampleList))) + "M"
    #     pool.apply_async(singleRun1, (refParams, tmpSample, mecatParams, optionTools, lncRNAcpatParams, dirSpec))
    # pool.close()
    # pool.join()
    # pybedtools.cleanup(remove_all=True)
    #
    # # singleRun1(refParams, allSampleMerged, mecatParams, optionTools, lncRNAcpatParams, dirSpec)
    #
    # for projectName in sampleDict:
    #     for tmpRun1, tmpRun2 in itertools.combinations(sampleDict[projectName], 2):
    #         from strainSpecificAS import getStrainSpecificComparion, specificAS
    #         getStrainSpecificComparion(tmpRun1, tmpRun2, dirSpec)

    # ngsCondFile = researchCfg.ngs_cfg_file
    # ngsCondDict = getSampleCond(ngsCondFile)
    for projectName in myproject2refStrain:
        for refStrain in myproject2refStrain[projectName]:
            # refParams = refInfoParams[refStrain]
            samples = myproject2refStrain[projectName][refStrain]
            # sampleNameList = [i.sampleName for i in samples]
            mergeProjects = []
            for s in samples:
                if "_".join(s.sampleName.split("_")[0:-1]) == refStrain:
                    mergeProjects.append(s)
            # from quantHybrid import quantifyHybrid
            # geneFeatureCountFile, transFeatureCountFile, sample2bamDict = quantifyHybrid(refParams=refParams, projects=samples, dirSpec=dirSpec, projectName=projectName, refStrain=refStrain)

            # from identifyDEhybrid import identifyDEhybrid2, deAnalysis2
            # compCondsDict = ngsCondDict[projectName][refStrain]
            # identifyDEhybrid(geneFeatureCountFile, ngsCondFile, refStrain, projectName, dirSpec,
            #                   compCondsDict=compCondsDict, sampleNameList=sampleNameList, analysisType="genes")
            # identifyDEhybrid(transFeatureCountFile, ngsCondFile, refStrain, projectName, dirSpec,
            #                   compCondsDict=compCondsDict, sampleNameList=sampleNameList, analysisType="transcripts")

            # identifyDEhybrid2(geneFeatureCountFile, sample2bamDict, samples, projectName, dirSpec, analysisType="gene")
            # identifyDEhybrid2(transFeatureCountFile, sample2bamDict, samples, projectName, dirSpec, analysisType="transcript")

            from identifyDAShybrid import identifyDAShybrid1, dasAnalysis, mergeIsoforms
            identifyDAShybrid1(mergeProjects, projectName, dirSpec)


def main():
    print("test start")
    defaultCfg = Config1(args.default_cfg)
    # integratedRun(defaultCfg=defaultCfg)
#    integratedRun1(defaultCfg=defaultCfg)
    test_integratedRun1(defaultCfg=defaultCfg)
    print("test end")

if __name__ == '__main__':
    main()
