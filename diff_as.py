#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: diff_as.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:20:35
Last modified: 2021-04-29 16:20:35
'''
from commonFuncs import *
from commonObjs import *
from itertools import combinations
import pandas as pd
import copy

def validateSamples(compCondFile, dataToProcess):
    samples = []
    conditions = set()
    if compCondFile:
        if validateFile(compCondFile):
            compCond = pd.read_csv(compCondFile, sep="\t")
            if "condition" in compCond.columns:
                for i, row in compCond.iterrows():
                    conditions.add(row.condition)
                    samples.extend([x for x in dataToProcess if x.project_name == row.project and x.condition == row.condition])
            else:
                raise Exception("Please set the condition in {}".format(compCondFile))
            samples = list(set(samples))
            # with open(compCondFile) as f:
            #     for line in f.readlines()[1:]:
            #         infoList = line.strip("\n").split("\t")
            #         samples = [x for x in dataToProcess if x.project_name == infoList[0] and x.condition == infoList[2]]
            #         conditions.add(infoList[2])
        else:
            compCondList = compCondFile.split(",")
            samples = [x for x in dataToProcess if x.condition in set(compCondList)]
            conditions = set(compCondList)
    else:
        samples = dataToProcess
    return samples, conditions

def mergeIsoforms(samples=None, dirSpec=None):
    tmpDict = {}
    mergedIso2ReadsBed = open("all_sample_merged_iso.bed", "w")
    for i in samples:
        isoformGroupedBed12 = os.path.join(dirSpec.out_dir, i.project_name, i.sample_name, "collapse", "isoformGrouped.bed12+")
        aseDir = os.path.join(dirSpec.out_dir, i.project_name, i.sample_name, "as_events", "ordinary_as")
        cmd = '''(cut -f 8,10 --output-delimiter=',' {}/PB/A3SS.confident.bed6+ {}/PB/A5SS.confident.bed6+ {}/PB/IR.confident.bed6+;
              cut -f 16,18 --output-delimiter=',' {}/PB/SE.confident.bed12+) | tr ',' '\n' | sort -u |
              filter.pl -o - {} -2 4 -m i > as_isoform.bed12+'''.format(aseDir, aseDir, aseDir, aseDir, isoformGroupedBed12)
        subprocess.call(cmd, shell=True)
        isoBedObj = BedFile("as_isoform.bed12+", type="bed12+")
        gene2iso = {}

        for iso in isoBedObj.reads:
            isoName = "{}_{}".format(i.sample_name, iso)
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

def getSample2Bam(samples, compCondFile, dirSpec=None):
    sample2bamDict = {"cond2bam": {}}
    compCondHeader = ["project", "sampleName", "condition", "repeat", "bamFile", "paired", "readsLength"]
    if os.path.isfile(compCondFile):
        compCond = pd.read_csv(compCondFile, sep="\t")
        if len(compCond.columns) == 7 and len(set(compCond.columns) & set(compCondHeader)) == 7:
            for i, row in compCond.iterrows():
                bamFile = row.bamFile
                if row.condition not in sample2bamDict["cond2bam"]:
                    sample2bamDict["cond2bam"][row.condition] = [(row.project + "_" + row.condition, bamFile, row.paired, row.readsLength)]
                else:
                    sample2bamDict["cond2bam"][row.condition].append((row.project + "_" + row.condition, bamFile, row.paired, row.readsLength))
        else:
            for dataObj in samples:
                alignDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "RNA-seq", "alignment")
                if dataObj.ngs_reads_paired == "paired":
                    leftReadsRepeats = dataObj.ngs_left_reads.split(";")
                    for i in range(len(leftReadsRepeats)):
                        repeatName = "repeat" + str(i)
                        bamFile = os.path.join(alignDir, repeatName, "{}.sorted.bam".format(repeatName))
                        if dataObj.condition not in sample2bamDict["cond2bam"]:
                            sample2bamDict["cond2bam"][dataObj.condition] = [(dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length)]
                        else:
                            sample2bamDict["cond2bam"][dataObj.condition].append((dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length))
                else:
                    if dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
                        singleReadsRepeats = dataObj.ngs_left_reads.split(";")
                    elif dataObj.ngs_left_reads == None and dataObj.ngs_right_reads:
                        singleReadsRepeats = dataObj.ngs_right_reads.split(";")
                    else:
                        raise Exception("The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

                    for i in range(len(singleReadsRepeats)):
                        repeatName = "repeat" + str(i)
                        bamFile = os.path.join(alignDir, repeatName, "{}.sorted.bam".format(repeatName))
                        if dataObj.condition not in sample2bamDict["cond2bam"]:
                            sample2bamDict["cond2bam"][dataObj.condition] = [(dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length)]
                        else:
                            sample2bamDict["cond2bam"][dataObj.condition].append((dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length))
    else:
        for dataObj in samples:
            alignDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "RNA-seq", "alignment")
            if dataObj.ngs_reads_paired == "paired":
                leftReadsRepeats = dataObj.ngs_left_reads.split(";")
                for i in range(len(leftReadsRepeats)):
                    repeatName = "repeat" + str(i)
                    bamFile = os.path.join(alignDir, repeatName, "{}.sorted.bam".format(repeatName))
                    if dataObj.condition not in sample2bamDict["cond2bam"]:
                        sample2bamDict["cond2bam"][dataObj.condition] = [(dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length)]
                    else:
                        sample2bamDict["cond2bam"][dataObj.condition].append((dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length))
            else:
                if dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
                    singleReadsRepeats = dataObj.ngs_left_reads.split(";")
                elif dataObj.ngs_left_reads == None and dataObj.ngs_right_reads:
                    singleReadsRepeats = dataObj.ngs_right_reads.split(";")
                else:
                    raise Exception("The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

                for i in range(len(singleReadsRepeats)):
                    repeatName = "repeat" + str(i)
                    bamFile = os.path.join(alignDir, repeatName, "{}.sorted.bam".format(repeatName))
                    if dataObj.condition not in sample2bamDict["cond2bam"]:
                        sample2bamDict["cond2bam"][dataObj.condition] = [(dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length)]
                    else:
                        sample2bamDict["cond2bam"][dataObj.condition].append((dataObj.sample_name + "_" + repeatName, bamFile, dataObj.ngs_reads_paired, dataObj.ngs_reads_length))
    return sample2bamDict

def diff_as1(dataToProcess, compCondFile=None, dirSpec=None, sampleMerged=False):
    print getCurrentTime() + " Identify differential alternative spliced isoforms..."
    if validateFile(compCondFile):
        compCondFile = os.path.abspath(compCondFile)
    else:
        compCondFile = compCondFile.split(",")
    prevDir = os.getcwd()
    dasDir = os.path.join(dirSpec.out_dir, "das")
    resolveDir(dasDir)
    samples, conditions = validateSamples(compCondFile, dataToProcess)
    if len(samples) == 0:
        raise Exception("May be you provide wrong sample comparison condition, please check it!")
    mergedIsoBed = mergeIsoforms(samples=samples, dirSpec=dirSpec)
    cmd = "bed2gpe.pl -g 13 {} > all_sample_merged_iso.gpe".format(mergedIsoBed)
    subprocess.call(cmd, shell=True)
    cmd = "genePredToGtf file all_sample_merged_iso.gpe all_sample_merged_iso.gtf"
    subprocess.call(cmd, shell=True)
    gtfFile = os.path.join(os.getcwd(), "all_sample_merged_iso.gtf")

    sample2bamDict = getSample2Bam(samples, compCondFile, dirSpec=dirSpec)


    currentThreads = 40
    for cond1, cond2 in combinations(conditions, 2):
        comp1Name = list(set([x[0] for x in sample2bamDict["cond2bam"][cond1]]))[0]
        comp2Name = list(set([x[0] for x in sample2bamDict["cond2bam"][cond2]]))[0]
        b1File = "b1File.txt"
        b2File = "b2File.txt"
        b1List = list(set([x[1] for x in sample2bamDict["cond2bam"][cond1]]))
        b2List = list(set([x[1] for x in sample2bamDict["cond2bam"][cond2]]))
        out1 = open(b1File, "w")
        print >> out1, ",".join(b1List)
        out1.close()
        out2 = open(b2File, "w")
        print >> out2, ",".join(b2List)
        out2.close()

        compOutDir = comp1Name + "_vs_" + comp2Name
        cond1Paired = list(set([x[2] for x in sample2bamDict["cond2bam"][cond1]]))
        cond2Paired = list(set([x[2] for x in sample2bamDict["cond2bam"][cond2]]))
        cond1Length = list(set([x[3] for x in sample2bamDict["cond2bam"][cond1]]))
        cond2Length = list(set([x[3] for x in sample2bamDict["cond2bam"][cond2]]))
        if cond1Paired != cond2Paired or len(cond1Paired) > 1 or len(cond2Paired) > 1 or cond1Length != cond2Length or len(cond1Length) > 1 or len(cond2Length) > 1:
            print "rMATS can't resolve the situation where condition {} and {} with different ngs sequencing strategy or read length!".format(comp1Name, comp2Name)
            continue
        cmd = "rmats.py --b1 {} --b2 {} --gtf {} --od {} -t {} --readLength {} --tstat {} --nthread {}"
        cmd = cmd.format(b1File, b2File, gtfFile, compOutDir, cond1Paired[0], cond1Length[0], currentThreads, currentThreads)
        subprocess.call(cmd, shell=True)

        resolveDir("{}.sigDiffAS".format(compOutDir), chdir=False)
        irDiff = pd.read_csv("{}/RI.MATS.JC.txt".format(compOutDir), sep="\t")
        seDiff = pd.read_csv("{}/SE.MATS.JC.txt".format(compOutDir), sep="\t")
        a3ssDiff = pd.read_csv("{}/A3SS.MATS.JC.txt".format(compOutDir), sep="\t")
        a5ssDiff = pd.read_csv("{}/A5SS.MATS.JC.txt".format(compOutDir), sep="\t")
        irDiff = irDiff.loc[irDiff.FDR <= 0.05]
        seDiff = seDiff.loc[seDiff.FDR <= 0.05]
        a3ssDiff = a3ssDiff.loc[a3ssDiff.FDR <= 0.05]
        a5ssDiff = a5ssDiff.loc[a5ssDiff.FDR <= 0.05]
        irDiff.to_csv("{}.sigDiffAS/IR.sig.txt".format(compOutDir), sep="\t", header=True, index=False)
        seDiff.to_csv("{}.sigDiffAS/SE.sig.txt".format(compOutDir), sep="\t", header=True, index=False)
        a3ssDiff.to_csv("{}.sigDiffAS/A3SS.sig.txt".format(compOutDir), sep="\t", header=True, index=False)
        a5ssDiff.to_csv("{}.sigDiffAS/A5SS.sig.txt".format(compOutDir), sep="\t", header=True, index=False)

    os.chdir(prevDir)
    print getCurrentTime() + " Identify differential alternative spliced isoforms done!"

# def diff_as(dataObj=None, refParams=None, dirSpec=None):
#     projectName, sampleName = dataObj.project_name, dataObj.sample_name
#     print getCurrentTime() + " Identify differential alternative spliced isoforms..."
#     prevDir = os.getcwd()
#     dasDir = os.path.join(dirSpec.out_dir, projectName, "hybrid", "das")
#     resolveDir(dasDir)
#
#     mergedIsoBed = mergeIsoforms(projects=samples, dirSpec=dirSpec)
#     cmd = "bed2gpe.pl -g 13 {} > all_sample_merged_iso.gpe".format(mergedIsoBed)
#     subprocess.call(cmd, shell=True)
#     cmd = "genePredToGtf file all_sample_merged_iso.gpe all_sample_merged_iso.gtf"
#     subprocess.call(cmd, shell=True)
#     gtfFile = os.path.join(os.getcwd(), "all_sample_merged_iso.gtf")
#
#     sample2bamDict = {"bamList": [], "bam2sample": {}, "sample2bam": {}}
#     for project in samples:
#         projectDir = os.path.join(dirSpec.out_dir, project.projectName, project.sampleName, "RNA-seq", "alignment")
#         if project.ngsReadPair == "paired":
#             leftReadsRepeats = project.ngsLeftReads.split(";")
#             for i in range(len(leftReadsRepeats)):
#                 repeatName = "repeat" + str(i)
#                 bamFile = os.path.join(projectDir, repeatName, "{}.sorted.bam".format(repeatName))
#                 sample2bamDict["bamList"].append(bamFile)
#                 sample2bamDict["bam2sample"][bamFile] = [project.sampleName + "_" + repeatName]
#                 sample2bamDict["sample2bam"][project.sampleName + "_" + repeatName] = [bamFile]
#         else:
#             if project.leftReads and project.rightReads == None:
#                 singleReadsRepeats = project.leftReads.split(";")
#             elif project.leftReads == None and project.rightReads:
#                 singleReadsRepeats = project.rightReads.split(";")
#             else:
#                 raise Exception(
#                     "The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")
#
#             for i in range(len(singleReadsRepeats)):
#                 repeatName = "repeat" + str(i)
#                 bamFile = os.path.join(projectDir, repeatName, "{}.sorted.bam".format(repeatName))
#                 sample2bamDict["bamList"].append(bamFile)
#                 sample2bamDict["bam2sample"][bamFile] = [project.sampleName + "_" + repeatName]
#                 sample2bamDict["sample2bam"][project.sampleName + "_" + repeatName] = [bamFile]
#
#     currentThreads = 40
#     for comp1, comp2 in combinations(samples, 2):
#         comp1Name = comp1.sampleName
#         comp2Name = comp2.sampleName
#         comp1repeats = len(comp1.ngsLeftReads.split(";"))
#         comp2repeats = len(comp2.ngsLeftReads.split(";"))
#         b1File = "b1File.txt"
#         b2File = "b2File.txt"
#         b1List = []
#         b2List = []
#         for i in range(comp1repeats):
#             repeatName = "repeat" + str(i)
#             b1List.extend(sample2bamDict["sample2bam"][comp1.sampleName + "_" + repeatName])
#         for i in range(comp2repeats):
#             repeatName = "repeat" + str(i)
#             b2List.extend(sample2bamDict["sample2bam"][comp2.sampleName + "_" + repeatName])
#         out1 = open(b1File, "w")
#         print >> out1, ",".join(b1List)
#         out1.close()
#         out2 = open(b2File, "w")
#         print >> out2, ",".join(b2List)
#         out2.close()
#
#         compOutDir = comp1Name + "_vs_" + comp2Name
#         if comp1.ngsPaired != comp2.ngsPaired or comp1.ngsReadsLength != comp2.ngsReadsLength >= 2:
#             print "rMATS can't resolve the situation where condition {} and {} in project {} with different ngs sequencing strategy or read length!".format(
#                 comp1Name, comp2Name, projectName)
#             continue
#         cmd = "rmats.py --b1 {} --b2 {} --gtf {} --od {} -t {} --readLength {} --tstat {} --nthread {}".format(b1File,
#                                                                                                                b2File,
#                                                                                                                gtfFile,
#                                                                                                                compOutDir,
#                                                                                                                comp1.ngsPaired,
#                                                                                                                comp1.ngsReadsLength,
#                                                                                                                currentThreads,
#                                                                                                                currentThreads)
#         subprocess.call(cmd, shell=True)
#
#     os.chdir(prevDir)
#     print getCurrentTime() + " Identify differential alternative spliced isoforms done!"
