#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: identifyDAShybrid.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-04 20:56:18
Last modified: 2019-12-04 20:56:19
'''

import argparse, datetime
from itertools import combinations
from multiprocessing import Pool
from Config import *
from commonObjs import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

##############################################
def dasAnalysis(comp, ngsCondDict, projectName, projects, dirSpec, hybridAlignDir, threads):
    cond1, cond2 = comp[0], comp[1]
    cond1SampleNames = ngsCondDict[cond2]
    cond2SampleNames = ngsCondDict[cond2]
    cond1Samples = [i for i in projects if i.sampleName in cond1SampleNames]
    cond2Samples = [i for i in projects if i.sampleName in cond2SampleNames]
    # group1 = ngsCondDict[projectName][cond1]["group"]
    # group2 = ngsCondDict[projectName][cond2]["group"]
    if len(cond1SampleNames) >= 2 and len(cond2SampleNames) >= 2:
        gtfFile = os.path.join(dirSpec.out_dir, projectName, "hybrid", "stringtie.merged.gtf")
        b1File = "b1File.txt"
        b2File = "b2File.txt"
        generateRMTAScompFile(hybridAlignDir, cond1SampleNames, b1File)
        generateRMTAScompFile(hybridAlignDir, cond2SampleNames, b2File)
        compOutDir = cond1 + "_vs_" + cond2
        try:
            group1ReadType = [i.ngsReadPair for i in cond1Samples]
            group2ReadType = [i.ngsReadPair for i in cond2Samples]
            group1ReadLength = [i.readLength for i in cond1Samples]
            group2ReadLength = [i.readLength for i in cond2Samples]
            if len(set(group1ReadType + group2ReadType)) >= 2 or len(set(group1ReadLength + group2ReadLength)) >= 2:
                print "rMATS can't resolve the situation where condition {} and {} in project {} with different ngs sequencing strategy or read length!".format(
                    cond1, cond2, projectName)
                return
            ngsReadPair = group1ReadType[0]
            ngsReadLength = group1ReadLength[0]
            cmd = "rmats.py --b1 {} --b2 {} --gtf {} --od {} -t {} --readLength {} --nthread {}".format(b1File, b2File, gtfFile, compOutDir, ngsReadPair, ngsReadLength, threads)
            subprocess.call(cmd, shell=True, executable="/bin/bash")
        except Exception as e:
            print e
    else:
        print "rMTAS can't analysis the comparison which any sample contains less than 2 repeats, so no differential alternative splicing analysis will be carried out!"
        # cmd = "EBseq.R {} {} {} {}".format(expFile, ngsCondFile, compGroups[i], compGroups[j])
        # subprocess.call(cmd, shell=True)

def identifyDAShybrid(projectName, projectsDict=None, ngsCondFile=None, refParams=None):
    print str(datetime.datetime.now()) + " Identify differential alternative spliced transcripts for project {}...".format(projectName)
    prevDir = os.getcwd()
    for projectName in projectsDict:
        hybridAlignDir = os.path.join(refParams.out_dir, projectName, "hybrid", "alignment")
        dasDir = os.path.join(refParams.out_dir, projectName, "hybrid", "das")
        resolveDir(dasDir)
        ngsCondDict = getSampleCond(ngsCondFile)
        if len(ngsCondDict[projectName]) == 1:
            print "Only one condition have been provided, so no differential alternative splicing analysis will be carried out!"
        else:
            compGroups = ngsCondDict[projectName]
            try:
                pool = Pool(processes=len(list(combinations(compGroups, 2))))
                multiResults = []
                for comp in list(combinations(compGroups, 2)):
                    singleRunRes = pool.apply_async(dasAnalysis, (comp, ngsCondDict, projectName, refParams, hybridAlignDir, projectsDict))
                    multiResults.append(singleRunRes)
                for j in multiResults:
                    j.wait()
            except Exception as e:
                print e
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Identify differential alternative spliced transcripts done!"

def mergeIsoforms(mergeProjects=None, projects=None, dirSpec=None):
    tmpDict = {}
    mergedIso2ReadsBed = open("all_sample_merged_iso.bed", "w")
    if not mergeProjects:
        mergeProjects = projects
    for i in mergeProjects:
        isoformGroupedBed12 = os.path.join(dirSpec.out_dir, i.projectName, i.sampleName, "postCorrEval", "isoformGrouped.bed12+")
        aseDir = os.path.join(dirSpec.out_dir, i.projectName, i.sampleName, "ASE")
        cmd = '''(cut -f 8,10 --output-delimiter=',' {}/PB/A3SS.confident.bed6+ {}/PB/A5SS.confident.bed6+ {}/PB/IR.bed6+;
              cut -f 16,18 --output-delimiter=',' {}/PB/SE.confident.bed12+) | tr ',' '\n' | sort -u |
              filter.pl -o - {} -2 4 -m i > as_isoform.bed12+'''.format(aseDir, aseDir, aseDir, aseDir, isoformGroupedBed12)
        subprocess.call(cmd, shell=True)
        isoBedObj = BedFile("as_isoform.bed12+", type="bed12+")
        gene2iso = {}

        for iso in isoBedObj.reads:
            isoName = "{}_{}".format(i.sampleName, iso)
            isoBedObj.reads[iso].name = isoName
            if isoBedObj.reads[iso].otherList[0] not in gene2iso:
                gene2iso[isoBedObj.reads[iso].otherList[0]] = []
            gene2iso[isoBedObj.reads[iso].otherList[0]].append(isoBedObj.reads[iso])

        for gene in gene2iso:
            for iso in gene2iso[gene]:
                if iso.chrom + "_" + iso.juncChain not in tmpDict:
                    tmpDict[iso.chrom + "_" + iso.juncChain] = [iso]
                else:
                    tmpDict[iso.chrom + "_" + iso.juncChain].append(iso)

    for tmp in tmpDict:
        isos = tmpDict[tmp]
        mergedIsoName = "+".join([x.name for x in isos])
        sortedIsos = sorted(isos, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
        repIso = copy.copy(sortedIsos[0])
        repIso.name = mergedIsoName
        print >> mergedIso2ReadsBed, str(repIso)
    mergedIso2ReadsBed.close()
    return os.path.join(os.getcwd(), "all_sample_merged_iso.bed")

def identifyDAShybrid1(samples, projectName, dirSpec):
    print str(datetime.datetime.now()) + " Identify differential alternative spliced transcripts for project {}...".format(projectName)
    prevDir = os.getcwd()
    dasDir = os.path.join(dirSpec.out_dir, projectName, "hybrid", "das")
    # gtfFile = os.path.join(dirSpec.out_dir, projectName, "hybrid", "stringtie.merged.gtf")
    resolveDir(dasDir)

    mergedIsoBed = mergeIsoforms(projects=samples, dirSpec=dirSpec)
    cmd = "bed2gpe.pl -g 13 {} > all_sample_merged_iso.gpe".format(mergedIsoBed)
    subprocess.call(cmd, shell=True)
    cmd = "genePredToGtf file all_sample_merged_iso.gpe all_sample_merged_iso.gtf"
    subprocess.call(cmd, shell=True)
    gtfFile = os.path.join(os.getcwd(), "all_sample_merged_iso.gtf")

    sample2bamDict = {"bamList": [], "bam2sample": {}, "sample2bam": {}}
    for project in samples:
        projectDir = os.path.join(dirSpec.out_dir, project.projectName, project.sampleName, "RNA-seq", "alignment")
        if project.ngsReadPair == "paired":
            leftReadsRepeats = project.ngsLeftReads.split(";")
            for i in range(len(leftReadsRepeats)):
                repeatName = "repeat" + str(i)
                bamFile = os.path.join(projectDir, repeatName, "{}.sorted.bam".format(repeatName))
                sample2bamDict["bamList"].append(bamFile)
                sample2bamDict["bam2sample"][bamFile] = [project.sampleName + "_" + repeatName]
                sample2bamDict["sample2bam"][project.sampleName + "_" + repeatName] = [bamFile]
        else:
            if project.leftReads and project.rightReads == None:
                singleReadsRepeats = project.leftReads.split(";")
            elif project.leftReads == None and project.rightReads:
                singleReadsRepeats = project.rightReads.split(";")
            else:
                raise Exception(
                    "The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

            for i in range(len(singleReadsRepeats)):
                repeatName = "repeat" + str(i)
                bamFile = os.path.join(projectDir, repeatName, "{}.sorted.bam".format(repeatName))
                sample2bamDict["bamList"].append(bamFile)
                sample2bamDict["bam2sample"][bamFile] = [project.sampleName + "_" + repeatName]
                sample2bamDict["sample2bam"][project.sampleName + "_" + repeatName] = [bamFile]

    currentThreads = 40
    for comp1, comp2 in combinations(samples, 2):
        comp1Name = comp1.sampleName
        comp2Name = comp2.sampleName
        comp1repeats = len(comp1.ngsLeftReads.split(";"))
        comp2repeats = len(comp2.ngsLeftReads.split(";"))
        b1File = "b1File.txt"
        b2File = "b2File.txt"
        b1List = []
        b2List = []
        for i in range(comp1repeats):
            repeatName = "repeat" + str(i)
            b1List.extend(sample2bamDict["sample2bam"][comp1.sampleName + "_" + repeatName])
        for i in range(comp2repeats):
            repeatName = "repeat" + str(i)
            b2List.extend(sample2bamDict["sample2bam"][comp2.sampleName + "_" + repeatName])
        out1 = open(b1File, "w")
        print >> out1, ",".join(b1List)
        out1.close()
        out2 = open(b2File, "w")
        print >> out2, ",".join(b2List)
        out2.close()

        compOutDir = comp1Name + "_vs_" + comp2Name
        if comp1.ngsPaired != comp2.ngsPaired or comp1.ngsReadsLength != comp2.ngsReadsLength >= 2:
            print "rMATS can't resolve the situation where condition {} and {} in project {} with different ngs sequencing strategy or read length!".format(
                comp1Name, comp2Name, projectName)
            continue
        cmd = "rmats.py --b1 {} --b2 {} --gtf {} --od {} -t {} --readLength {} --tstat {} --nthread {}".format(b1File, b2File, gtfFile, compOutDir, comp1.ngsPaired, comp1.ngsReadsLength, currentThreads, currentThreads)
        subprocess.call(cmd, shell=True)

    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Identify differential alternative spliced transcripts done!"

def identifyDAShybrid2(projectName, dirSpec, compCondsDict=None, projects=None):
    print str(datetime.datetime.now()) + " Identify differential alternative spliced transcripts for project {}...".format(projectName)
    prevDir = os.getcwd()
    hybridAlignDir = os.path.join(dirSpec.out_dir, projectName, "hybrid", "alignment")
    dasDir = os.path.join(dirSpec.out_dir, projectName, "hybrid", "das")
    resolveDir(dasDir)
    conditions = compCondsDict.keys()
    currentThreads = sum([i.threads for i in projects]) if sum([i.threads for i in projects]) < 64 else 64
    try:
        pool = Pool(processes=len(list(combinations(conditions, 2))))
        for comp in list(combinations(conditions, 2)):
            pool.apply_async(dasAnalysis, (comp, compCondsDict, projectName, projects, dirSpec, hybridAlignDir, currentThreads))
        pool.close()
        pool.join()
    except Exception as e:
        print e
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Identify differential alternative spliced transcripts done!"

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    ngsCondFile = researchCfg.ngs_sample_cond
    identifyDAShybrid(projectsDict=projectsDict, ngsCondFile=ngsCondFile, refParams=refParams)

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)
