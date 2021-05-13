#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: collapse.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:16:07
Last modified: 2021-04-29 16:16:07
'''
from commonFuncs import *
from commonObjs import *
import pybedtools

def readsAssign(transBedFile, readsBedFile, offset=10, minConsenesusIntronN=1, minCoverageOnRead=0.9, singleLine=True, transColNum=13, readsColNum=13, outPrefix="readsAssign", group=True):
    transBedObj = pybedtools.BedTool(transBedFile)
    readsBedObj = pybedtools.BedTool(readsBedFile)
    intersectRes = readsBedObj.intersect(transBedObj, wa=True, wb=True)
    notintersectRes = readsBedObj.intersect(transBedObj, v=True)
    matchedDict = {}
    outList = []
    for i in intersectRes:
        infoList = str(i).strip("\n").split("\t")
        if readsColNum > 12:
            readsBed = Bed12Plus("\t".join(infoList[:readsColNum]))
        else:
            readsBed = Bed12("\t".join(infoList[:readsColNum]))
        if transColNum < 13:
            raise Exception("The colmun count of your transcript bed should large than 12")
        transBed = Bed12Plus("\t".join(infoList[readsColNum:readsColNum+transColNum]))
        readStart, readEnd, readsExons, transExons = readsBed.chromStart, readsBed.chromEnd, readsBed.exons, transBed.exons
        overlapLength = getBlockLength(getOverlapOfTuple(readsExons, transExons))
        coverageOnRead = overlapLength/float(getBlockLength(readsExons))
        coverageOnTrans = overlapLength/float(getBlockLength(transExons))
        transGene = transBed.record[12]
        if len(readsBed.exons) == 1:
            if transBed.exonStarts[0] < readEnd and readEnd <= transBed.exonEnds[0] or \
                    transBed.exonStarts[-1] <= readStart and readStart < transBed.exonEnds[-1]:
                inEndExonAndExtension = 1
            else:
                inEndExonAndExtension = 0

            if readsBed.name not in matchedDict:
                matchedDict[readsBed.name] = {"matchedTrans": [transBed.name], "matchedCovOnReads": [coverageOnRead],
                                              "matchedCovOnTrans": [coverageOnTrans], "matchedGenes": [transGene],
                                              "matchedTransExonCount": [len(transBed.exons)],
                                              "inEndExonAndExtension": [inEndExonAndExtension], "readsBed": readsBed}
            else:
                matchedDict[readsBed.name]["matchedTrans"].append(transBed.name)
                matchedDict[readsBed.name]["matchedCovOnReads"].append(coverageOnRead)
                matchedDict[readsBed.name]["matchedCovOnTrans"].append(coverageOnTrans)
                matchedDict[readsBed.name]["matchedTransExonCount"].append(len(transBed.exons))
                matchedDict[readsBed.name]["matchedGenes"].append(transGene)
                matchedDict[readsBed.name]["inEndExonAndExtension"].append(inEndExonAndExtension)
        else:
            consensusIntronN = getConsensusIntronN(readsBed.exonStarts, readsBed.exonEnds, transBed.exonStarts, transBed.exonEnds, offset=offset)
            readsJuncChain, transJuncChain = readsBed.juncChain, transBed.juncChain
            if re.search(readsJuncChain + "$", transJuncChain) or re.search("^" + readsJuncChain, transJuncChain):
                juncChainFlank = 1
            else:
                juncChainFlank = 0

            if readsBed.name not in matchedDict:
                matchedDict[readsBed.name] = {"matchedTrans": [transBed.name], "matchedCovOnReads": [coverageOnRead],
                                              "matchedCovOnTrans": [coverageOnTrans], "matchedGenes": [transGene],
                                              "consensusIntronN": [consensusIntronN],
                                              "juncChainFlank": [juncChainFlank], "readsBed": readsBed}
            else:
                matchedDict[readsBed.name]["matchedTrans"].append(transBed.name)
                matchedDict[readsBed.name]["matchedCovOnReads"].append(coverageOnRead)
                matchedDict[readsBed.name]["matchedCovOnTrans"].append(coverageOnTrans)
                matchedDict[readsBed.name]["matchedGenes"].append(transGene)
                matchedDict[readsBed.name]["consensusIntronN"].append(consensusIntronN)
                matchedDict[readsBed.name]["juncChainFlank"].append(juncChainFlank)

    for readName in matchedDict:
        if singleLine:
            if "juncChainFlank" in matchedDict[readName]:
                readType = "I"
                readIntronN = len(matchedDict[readName]["readsBed"].introns)
                realMatchedTrans, realConsensusIntronN, realMatchedCovOnReads, realMatchedCovOnTrans, realMatchedGene = [], [], [], [], []
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    if matchedDict[readName]["juncChainFlank"][j] == 1 or \
                        (readIntronN < minConsenesusIntronN and readIntronN == matchedDict[readName]["consensusIntronN"][j]) or \
                        (readIntronN >= minConsenesusIntronN and matchedDict[readName]["consensusIntronN"][j] >= minConsenesusIntronN) or \
                        (matchedDict[readName]["matchedCovOnReads"][j] >= minCoverageOnRead):
                        readType = "E"
                        realMatchedTrans.append(matchedDict[readName]["matchedTrans"][j])
                        realConsensusIntronN.append(matchedDict[readName]["consensusIntronN"][j])
                        realMatchedCovOnReads.append(matchedDict[readName]["matchedCovOnReads"][j])
                        realMatchedCovOnTrans.append(matchedDict[readName]["matchedCovOnTrans"][j])
                        realMatchedGene.append(matchedDict[readName]["matchedGenes"][j])
                if readType == "E":
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(realMatchedTrans), ",".join(realMatchedGene),
                                    ",".join(map(str, realConsensusIntronN)), ",".join(map(str, realMatchedCovOnReads)),
                                    ",".join(map(str, realMatchedCovOnTrans))])
                else:
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(matchedDict[readName]["matchedTrans"]),
                                    ",".join(matchedDict[readName]["matchedGenes"]),
                                    ",".join(map(str, matchedDict[readName]["consensusIntronN"])),
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnReads"])),
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnTrans"]))])
            else:
                readType = "I"
                realMatchedTrans, realMatchedCovOnReads, realMatchedCovOnTrans, realMatchedGene = [], [], [], []
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    if matchedDict[readName]["matchedTransExonCount"][j] > 1 and \
                        matchedDict[readName]["inEndExonAndExtension"][j] == 0 and \
                        matchedDict[readName]["matchedCovOnReads"][j] < minCoverageOnRead:
                        continue
                    readType = "E"
                    realMatchedTrans.append(matchedDict[readName]["matchedTrans"][j])
                    realMatchedCovOnReads.append(matchedDict[readName]["matchedCovOnReads"][j])
                    realMatchedCovOnTrans.append(matchedDict[readName]["matchedCovOnTrans"][j])
                    realMatchedGene.append(matchedDict[readName]["matchedGenes"][j])
                if readType == "E":
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(map(str, realMatchedTrans)), ",".join(map(str, realMatchedGene)), "NA",
                                    ",".join(map(str, realMatchedCovOnReads)), ",".join(map(str, realMatchedCovOnTrans))])
                else:
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(matchedDict[readName]["matchedTrans"]),
                                    ",".join(matchedDict[readName]["matchedGenes"]), "NA",
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnReads"])),
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnTrans"]))])
        else:
            if "juncChainFlank" in matchedDict[readName]:
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    readType = "I"
                    exonNum = len(matchedDict[readName]["readsBed"].exons)
                    if matchedDict[readName]["juncChainFlank"][j] == 1 or \
                            (exonNum < minConsenesusIntronN and exonNum - 1 == matchedDict[readName]["consensusIntronN"][j]) or \
                            (exonNum >= minConsenesusIntronN and matchedDict[readName]["consensusIntronN"][j] >= minConsenesusIntronN):
                        readType = "E"
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    matchedDict[readName]["matchedTrans"][j], matchedDict[readName]["matchedGenes"][j],
                                    matchedDict[readName]["consensusIntronN"][j],
                                    matchedDict[readName]["matchedCovOnReads"][j],
                                    matchedDict[readName]["matchedCovOnTrans"][j]])
            else:
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    if matchedDict[readName]["matchedTransExonCount"][j] > 1 and \
                        matchedDict[readName]["inEndExonAndExtension"][j] == 0 and \
                        matchedDict[readName]["matchedCovOnReads"][j] < minCoverageOnRead:
                        readType = "I"
                    else:
                        readType = "E"
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    matchedDict[readName]["matchedTrans"][j],
                                    matchedDict[readName]["matchedCovOnReads"][j], "NA",
                                    matchedDict[readName]["matchedTrans"][j],
                                    matchedDict[readName]["matchedCovOnTrans"][j]])

    for i in notintersectRes:
        outList.append([str(i).strip("\n"), "IG", "NA", "NA", "NA", "NA", "NA"])
        # print str(i).strip("\n") + "\t" + "\t".join(["IG", "NA", "NA", "NA", "NA", "NA"])

    assignOut = open(outPrefix + ".bed12+", "w")
    sortedOutList = sorted(outList, key=lambda x: (x[0].split("\t")[0], int(x[0].split("\t")[1])))
    for item in sortedOutList:
        print >> assignOut, "\t".join(map(str, item))
    assignOut.close()

    if group:
        unambiOut = open(outPrefix + ".unambi.bed12+", "w")
        novelList = []
        for item in sortedOutList:
            if item[3] == "NA" or item[1] != "E":
                novelList.append(item[0].split("\t"))
            else:
                uniqGenes = list(set(item[3].split(",")))
                if len(uniqGenes) == 1:
                    print >> unambiOut, item[0] + "\t" + uniqGenes[0]
                else:
                    novelList.append(item[0].split("\t"))
        sortedNovelList = sorted(novelList, key=lambda x: int(x[1]))
        inc, clusterEnd = 0, 0
        for novelItem in sortedNovelList:
            myStart, myEnd = int(novelItem[1]), int(novelItem[2])
            if myStart < clusterEnd:
                if myEnd > clusterEnd:
                    clusterEnd = myEnd
            else:
                inc += 1
                clusterEnd = myEnd
            print >> unambiOut, "\t".join(novelItem) + "\t" + ":".join(map(str, [novelItem[0], novelItem[5], inc]))
        unambiOut.close()

def collapse(dataObj=None, collapseParams=None, refParams=None, dirSpec=None, threads=10):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Collapse with cDNA_cupcake for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)

    resolveDir(baseDir, "collapse")
    logDir = os.path.join(baseDir, "log")
    resolveDir(logDir, chdir=False)
    processedFa = os.path.join(baseDir, "filtration", "processed.fa")
    processedBed = os.path.join(baseDir, "filtration", "processed.bed12+")
    flncSam = os.path.join(baseDir, "mapping", "flnc.mm2.sam")
    processedIds = getFxSequenceId(processedFa, isFa=True)
    getSubSamByName(flncSam, nameList=processedIds, isBam=False, nameListIsFile=False, outPrefix="processed", sort=True,
                    threads=threads)
    if collapseParams.dun_merge_5_shorter:
        cmd = "collapse_isoforms_by_sam.py --input {} -s processed.sorted.sam --max_5_diff {} --max_3_diff {} " \
              "--flnc_coverage {} -c {} -i {} --max_fuzzy_junction {} --dun-merge-5-shorter -o tofu 1>{}/tofu.collapse.log 2>&1"
    else:
        cmd = "collapse_isoforms_by_sam.py --input {} -s processed.sorted.sam --max_5_diff {} --max_3_diff {} " \
              "--flnc_coverage {} -c {} -i {} --max_fuzzy_junction {} -o tofu 1>{}/tofu.collapse.log 2>&1"
    cmd = cmd.format(processedFa, collapseParams.max_5_diff, collapseParams.max_3_diff, collapseParams.fl_coverage,
                     collapseParams.coverage, collapseParams.identity, collapseParams.max_fuzzy_junction, logDir)
    subprocess.call(cmd, shell=True)

    cmd = "gtfToGenePred tofu.collapsed.good.gff tofu.collapsed.gpe -genePredExt"
    subprocess.call(cmd, shell=True)
    cmd = "gpe2bed.pl tofu.collapsed.gpe -g > tofu.collapsed.bed12+"
    subprocess.call(cmd, shell=True)
    readsAssign(refParams.ref_bed, "tofu.collapsed.bed12+", readsColNum=13, outPrefix="tofu.collapsed.assigned",
                group=True)
    cmd = "cut -f1-12,14 tofu.collapsed.assigned.unambi.bed12+ > isoformGrouped.bed12+"
    subprocess.call(cmd, shell=True)

    cmd = '''seqkit grep {} -f <(cut -f 2 tofu.collapsed.group.txt | tr ',' '\n') -w 0 > processed.ignore_id_removed.fa'''.format(processedFa)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''filter.pl -o <(cut -f 2 tofu.collapsed.group.txt | tr ',' '\n') {} -2 4 -m i > processed.ignore_id_removed.bed12+'''.format(processedBed)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    readsAssign(refParams.ref_bed, "processed.ignore_id_removed.bed12+", readsColNum=14, outPrefix="reads.assigned", group=True)

    cmd = "cut -f 1-12,15 reads.assigned.unambi.bed12+ | bed2gpe.pl -b 12 -g 13 - | genePredToGtf file stdin reads.unambi.gtf -source=iFLAS"
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print getCurrentTime() + " Collapse with cDNA_cupcake for project {} sample {} done!".format(projectName, sampleName)
