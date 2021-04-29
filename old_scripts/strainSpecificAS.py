#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: strainSpecificAS.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-05-17 12:57:09
Last modified: 2020-05-17 12:57:09
'''

import argparse
from commonObjs import *
from Config import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

##############################################
def processEsFile(esFile):
    myDict = {}
    recordToExon = {}
    with open(esFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            chrom = records[0]
            esInfoList = records[3].split("@")
            esExons = esInfoList[1].split(";")
            recordToExon["\t".join(records)] = []
            for i in esExons:
                exon = "{}:{}".format(chrom, i)
                if exon not in myDict:
                    myDict[exon] = [records]
                else:
                    myDict[exon].append(records)
                recordToExon["\t".join(records)].append(exon)
    return myDict, recordToExon

def getAsEvents(bedFile):
    myDict = {}
    for r in BedFile(bedFile, type="bed6").reads:
        asEvent = r.chrom + ":" + r.strand + ":" + r.name.split(":")[-1]
        myDict[asEvent] = r.record
    return myDict

def processBedFile(bedFile):
    myDict = {}
    with open(bedFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            chrom, chromStart = records[0], int(records[1])
            blockSizes = [int(i) for i in records[10].split(",")]
            exonStarts = [chromStart + int(i) for i in records[11].split(",")]
            exonEnds = [exonStarts[i] + blockSizes[i] for i in range(len(blockSizes))]
            for i in range(len(exonStarts)):
                exon = "{}:{}-{}".format(chrom, exonStarts[i], exonEnds[i])
                if exon not in myDict:
                    myDict[exon] = ""
    return myDict

def filter(fref=None, source=None, target=None, source_all=None, target_all=None, mode=None, outFile=None):
    out = open(outFile, "w")
    refExonDict = {}
    with open(fref) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            transName, chrom = records[0], records[1]
            exonStarts = [int(i) for i in records[8].strip(',').split(',')]
            exonEnds = [int(i) for i in records[9].strip(',').split(',')]
            for i in range(len(exonStarts)):
                exon = "{}:{}-{}".format(chrom, exonStarts[i], exonEnds[i])
                if exon not in refExonDict:
                    refExonDict[exon] = [transName]
                else:
                    refExonDict[exon].append(transName)
    sourceDict, sourceRecordToExonDict = processEsFile(source)
    targetDict, targetRecordToExonDict = processEsFile(target)
    sourceAllDict = processBedFile(source_all)
    targetAllDict = processBedFile(target_all)
    refExons = set(refExonDict.keys())
    sourceExons = set(sourceDict.keys())
    sourceAllExons = set(sourceAllDict.keys())
    targetExons = set(targetDict.keys())
    targetAllExons = set(targetAllDict.keys())
    if mode == "e1":
        filterExons = targetExons - (sourceAllExons | refExons)
        for i in filterExons:
            for j in targetDict[i]:
                if len(set(targetRecordToExonDict["\t".join(j)]) & (sourceAllExons | refExons)) == 0:
                    print >> out, "\t".join(j)
    elif mode == "e2":
        filterExons = sourceExons - targetAllExons
        for i in filterExons:
            if i in sourceDict:
                for j in sourceDict[i]:
                    if len(set(sourceRecordToExonDict["\t".join(j)]) & targetAllExons) == 0:
                        print >> out, "\t".join(j)
    else:
        filterExons = targetExons & (refExons | sourceExons)
        for i in filterExons:
            for j in targetDict[i]:
                print >>out, "\t".join(j)
    out.close()

def specificAS(comp1AsFile, comp2AsFile, asType="SE", comp1RefAsFile=None, comp2RefAsFile=None, mode="", offset=0, outFile=None):
    print str(datetime.datetime.now()) + " Get strain/species-specific alternative splice events..."

    comp1AsEvents = getAsEvents(comp1AsFile)
    comp2AsEvents = getAsEvents(comp2AsFile)
    # isoform1Dict = BedFile(isoform1File, type="bed12")
    # isoform2Dict = BedFile(isoform2File, type="bed12")

    if comp1RefAsFile:
        comp1RefAsEvents = getAsEvents(comp1RefAsFile)
        comp1AsEvents.update(comp1RefAsEvents)
    if comp2RefAsFile:
        comp2RefAsEvents = getAsEvents(comp2RefAsFile)
        comp2AsEvents.update(comp2RefAsEvents)

    commonEvents = []
    if asType == "SE":
        for event1 in comp1AsEvents:
            event1Split = event1.split(":")
            chrom1, strand1 = event1Split[0:2]
            left1, right1 = int(event1Split[2].split("@")[0]), int(event1Split[2].split("@")[-1])
            skipedExon1 = event1Split[2].split("@")[1]
            for event2 in comp2AsEvents:
                event2Split = event2.split(":")
                chrom2, strand2 = event2Split[0:2]
                left2, right2 = int(event2Split[2].split("@")[0]), int(event2Split[2].split("@")[-1])
                skipedExon2 = event2Split[2].split("@")[1]
                if chrom1 == chrom2 and strand1 == strand2 and skipedExon1 == skipedExon2:
                    if abs(left1 - left2) <= offset and abs(right1 - right2) <= offset:
                        if mode == "e1":
                            commonEvents.append(event1)
                        else:
                            commonEvents.append(event2)

    else:
        for event1 in comp1AsEvents:
            event1Split = event1.split(":")
            chrom1, strand1 = event1Split[0:2]
            left1, right1 = listStr2Int(event1Split[2].split("-"))
            for event2 in comp2AsEvents:
                event2Split = event2.split(":")
                chrom2, strand2 = event2Split[0:2]
                left2, right2 = listStr2Int(event2Split[2].split("-"))
                if chrom1 == chrom2 and strand1 == strand2:
                    if abs(left1 - left2) <= offset and abs(right1 - right2) <= offset:
                        if mode == "e1":
                            commonEvents.append(event1)
                        else:
                            commonEvents.append(event2)

    if mode == "e1":
        specificEvents = set(comp1AsEvents.keys()) - set(commonEvents)
        for i in specificEvents:
            print >>outFile, "\t".join(comp1AsEvents[i])
    else:
        specificEvents = set(comp2AsEvents.keys()) - set(commonEvents)
        for i in specificEvents:
            print >>outFile, "\t".join(comp2AsEvents[i])

    print str(datetime.datetime.now()) + " Get strain/species-specific alternative splice events done!"

def getStrainSpecificComparion(comp1, comp2, dirSpec):
    tgsAnalysisSet = set([])
    asTypes = ["SE", "IR", "A5SS", "A3SS"]
    if comp1.projectName == comp2.projectName:
        strainSpecificAnalysisDir = os.path.join(dirSpec.out_dir, comp1.projectName, "strainSpecific")
        prevDir = os.getcwd()
        resolveDir(strainSpecificAnalysisDir)
        comp1dir = os.path.join(dirSpec.out_dir, comp1.projectName, comp1.sampleName)
        comp2dir = os.path.join(dirSpec.out_dir, comp2.projectName, comp1.sampleName)

        compName = comp1.sampleName + "vs" + comp2.sampleName

        # isoform1File = os.path.join(comp1dir, "ASE", "deSingleExonRead.AS.confident.bed12+")
        # isoform2File = os.path.join(comp2dir, "ASE", "deSingleExonRead.AS.confident.bed12+")

        for asType in asTypes:
            if asType == "SE":
                asFile = asType + ".bed12+"
                refAsFile = asType + ".reference.bed12+"
            else:
                asFile = asType + ".bed6+"
                refAsFile = asType + ".reference.bed6+"
            comp1AsFile = os.path.join(comp1dir, "ASE", "PB", asFile)
            comp2AsFile = os.path.join(comp2dir, "ASE", "PB", asFile)
            comp1RefAsFile = os.path.join(comp1dir, "ASE", "characterization", refAsFile)
            comp2RefAsFile = os.path.join(comp2dir, "ASE", "characterization", refAsFile)

            comp1SpecificOutFile = ".".join([compName, comp1.sampleName, "specific", asType, "txt"])
            comp2SpecificOutFile = ".".join([compName, comp2.sampleName, "specific", asType, "txt"])
            specificAS(comp1AsFile, comp2AsFile, asType=asType,
                       comp1RefAsFile=comp1RefAsFile, comp2RefAsFile=comp2RefAsFile, mode="e1", offset=3,
                       outFile=comp1SpecificOutFile)
            specificAS(comp1AsFile, comp2AsFile, asType=asType,
                       comp1RefAsFile=comp1RefAsFile, comp2RefAsFile=comp2RefAsFile, mode="e2", offset=3,
                       outFile=comp2SpecificOutFile)
        os.chdir(prevDir)

    # for item in itertools.combinations(projects, 2):
    #     if item[0].group != item[1].group and item[0].projectName == item[1].projectName:
    #         strainSpecificAnalysisDir = os.path.join(dirSpec.out_dir, item[0].projectName, "strainSpecific")
    #         prevDir = os.getcwd()
    #         resolveDir(strainSpecificAnalysisDir)
    #         if item[0].tgsDataDir not in tgsAnalysisSet or item[1].tgsDataDir not in tgsAnalysisSet:
    #             tgsAnalysisSet.add(item[0].tgsDataDir)
    #             tgsAnalysisSet.add(item[1].tgsDataDir)
    #
    #             strain1dir = os.path.join(dirSpec.out_dir, item[0].projectName, item[0].group)
    #             strain2dir = os.path.join(dirSpec.out_dir, item[1].projectName, item[1].group)
    #             isoform1File = os.path.join(strain1dir, "ASE", "deSingleExonRead.AS.confident.bed12+")
    #             isoform2File = os.path.join(strain2dir, "ASE", "deSingleExonRead.AS.confident.bed12+")
    #
    #             for asType in asTypes:
    #                 if asType == "SE":
    #                     asFile = asType + ".bed12+"
    #                     refAsFile = asType + ".reference.bed12+"
    #                 else:
    #                     asFile = asType + ".bed6+"
    #                     refAsFile = asType + ".reference.bed6+"
    #                 strain1AsFile = os.path.join(strain1dir, "ASE", "PB", asFile)
    #                 strain2AsFile = os.path.join(strain2dir, "ASE", "PB", asFile)
    #                 strain1RefAsFile = os.path.join(strain1dir, "ASE", "characterization", refAsFile)
    #                 strain2RefAsFile = os.path.join(strain2dir, "ASE", "characterization", refAsFile)
    #
    #                 strain1SpecificOutFile = ".".join([item[0].group, "specific", asType, "txt"])
    #                 strain2SpecificOutFile = ".".join([item[1].group, "specific", asType, "txt"])
    #                 strainSpecificAS(strain1AsFile, isoform1File, strain2AsFile, isoform2File, asType=asType,
    #                                  strain1RefAS=strain1RefAsFile, strain2RefAS=strain2RefAsFile, mode="e1", offset=3,
    #                                  outFile=strain1SpecificOutFile)
    #                 strainSpecificAS(strain1AsFile, isoform1File, strain2AsFile, isoform2File, asType=asType,
    #                                  strain1RefAS=strain1RefAsFile, strain2RefAS=strain2RefAsFile, mode="e2", offset=3,
    #                                  outFile=strain2SpecificOutFile)
    #         os.chdir(prevDir)

def getStrainSpecificComparion1(dirSpec, projects=None):
    asTypes = ["SE", "IR", "A5SS", "A3SS"]
    for item in itertools.combinations(projects, 2):
        alleleSpecificAnalysisDir = os.path.join(dirSpec.out_dir, item[0].projectName, "alleleSpecific")
        prevDir = os.getcwd()
        resolveDir(alleleSpecificAnalysisDir)
        # if item[0].tgsDataDir not in tgsAnalysisSet or item[1].tgsDataDir not in tgsAnalysisSet:
        #     tgsAnalysisSet.add(item[0].tgsDataDir)
        #     tgsAnalysisSet.add(item[1].tgsDataDir)

        strain1dir = os.path.join(dirSpec.out_dir, item[0].projectName, item[0].sampleName)
        strain2dir = os.path.join(dirSpec.out_dir, item[1].projectName, item[1].sampleName)
        isoform1File = os.path.join(strain1dir, "ASE", "deSingleExonRead.AS.confident.bed12+")
        isoform2File = os.path.join(strain2dir, "ASE", "deSingleExonRead.AS.confident.bed12+")

        for asType in asTypes:
            if asType == "SE":
                asFile = asType + ".bed12+"
                refAsFile = asType + ".reference.bed12+"
            else:
                asFile = asType + ".bed6+"
                refAsFile = asType + ".reference.bed6+"
            strain1AsFile = os.path.join(strain1dir, "ASE", "PB", asFile)
            strain2AsFile = os.path.join(strain2dir, "ASE", "PB", asFile)
            strain1RefAsFile = os.path.join(strain1dir, "ASE", "characterization", refAsFile)
            strain2RefAsFile = os.path.join(strain2dir, "ASE", "characterization", refAsFile)

            strain1SpecificOutFile = ".".join([item[0].sampleName, "specific", asType, "txt"])
            strain2SpecificOutFile = ".".join([item[1].sampleName, "specific", asType, "txt"])
            strainSpecificAS(strain1AsFile, isoform1File, strain2AsFile, isoform2File, asType=asType,
                             strain1RefAS=strain1RefAsFile, strain2RefAS=strain2RefAsFile, mode="e1", offset=3,
                             outFile=strain1SpecificOutFile)
            strainSpecificAS(strain1AsFile, isoform1File, strain2AsFile, isoform2File, asType=asType,
                             strain1RefAS=strain1RefAsFile, strain2RefAS=strain2RefAsFile, mode="e2", offset=3,
                             outFile=strain2SpecificOutFile)
        os.chdir(prevDir)

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    getStrainSpecificComparion(refParams=refParams, projects=projects)

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)
