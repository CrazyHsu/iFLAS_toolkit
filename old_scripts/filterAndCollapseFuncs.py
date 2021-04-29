#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: filterAndCollapseFuncs.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-08 17:31:43
Last modified: 2019-12-08 19:05:28
'''
import datetime, hashlib, pybedtools, pysam
from commonFuncs import *
from commonObjs import *

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
            consensusIntronN = getConsensusIntronN1(readsBed.exonStarts, readsBed.exonEnds, transBed.exonStarts, transBed.exonEnds, offset=offset)
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

    # if group:
    #     unambiOut = open(outPrefix + ".unambi.bed12+", "w")
    #     ambiOut = open(outPrefix + ".ambiguous.bed12+", "w")
    #     novelList = []
    #     for item in sortedOutList:
    #         if item[3] == "NA" or item[1] != "E":
    #             novelList.append(item[0].split("\t"))
    #         else:
    #             uniqGenes = list(set(item[3].split(",")))
    #             if len(uniqGenes) == 1:
    #                 print >> unambiOut, item[0] + "\t" + uniqGenes[0]
    #             else:
    #                 print >> ambiOut, item[0] + "\t" + ",".join(uniqGenes)
    #     sortedNovelList = sorted(novelList, key=lambda x: int(x[1]))
    #     inc, clusterEnd = 0, 0
    #     for novelItem in sortedNovelList:
    #         myStart, myEnd = int(novelItem[1]), int(novelItem[2])
    #         if myStart < clusterEnd:
    #             if myEnd > clusterEnd:
    #                 clusterEnd = myEnd
    #         else:
    #             inc += 1
    #             clusterEnd = myEnd
    #         print >> unambiOut, "\t".join(novelItem) + "\t" + ":".join(map(str, [novelItem[0], novelItem[5], inc]))
    #     unambiOut.close()
    #     ambiOut.close()

    if group:
        unambiOut = open(outPrefix + ".unambi.bed12+", "w")
        # ambiOut = open(outPrefix + ".ambiguous.bed12+", "w")
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
                    # print >> ambiOut, item[0] + "\t" + ",".join(uniqGenes)
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
        # ambiOut.close()

def strandAdjust(genomeFasta, refGPE, bedFile, minCoverage, juncDiffScore, strandAdjust=None, strandConfirmed=None):
    tmpBed = "tmp.bed"
    GenePredObj(refGPE).toBed(outFile=tmpBed)

    strandAdjustOut = open(strandAdjust, "w")
    strandConfirmedOut = open(strandConfirmed, "w")

    bedObj = BedFile(bedFile, type="bed12+")
    gpeBedObj = pybedtools.bedtool.BedTool(tmpBed)
    singleExonReadBedList = []
    # multiExonReadDinucleotideDict = {}
    multiExonReadDinucleotideBedList = []
    for readName in bedObj.reads:
        # multiExonReadDinucleotideDict[readName] = {}
        read = bedObj.reads[readName]
        chrom, start, end, strand = read.chrom, read.chromStart, read.chromEnd, read.strand
        readExonStarts, readExonEnds = read.exonStarts, read.exonEnds
        # readExons = [(readExonStarts[x], readExonEnds[x]) for x in range(len(readExonStarts))]
        readRegionBed = " ".join(map(str, [chrom, start, end, readName, ".", strand]))
        blockCount = read.blockCount
        if blockCount == 1:
            singleExonReadBedList.append(readRegionBed)
        else:
            readJuncStarts, readJuncEnds = getIntrons(readExonStarts, readExonEnds)
            for i in range(len(readJuncStarts)):
                leftDinucleotidePos = " ".join(map(str, [chrom, readJuncStarts[i], readJuncStarts[i] + 2, ":".join([readName, "junction"+str(i), "left"]), ".", strand]))
                rightDinucleotidePos = " ".join(map(str, [chrom, readJuncEnds[i] - 2, readJuncEnds[i], ":".join([readName, "junction"+str(i), "right"]), ".", strand]))
                multiExonReadDinucleotideBedList.extend([leftDinucleotidePos, rightDinucleotidePos])
                # multiExonReadDinucleotideDict[readName].update({"junction"+str(i): {"left": leftDinucleotidePos, "right": rightDinucleotidePos}})

    # single exon process
    singleExonReadBedObj = pybedtools.bedtool.BedTool("\n".join(singleExonReadBedList), from_string=True)
    singleExonIntersectRes = gpeBedObj.intersect(singleExonReadBedObj, wa=True, wb=True)
    singleExonStrandDict = {}
    for i in singleExonIntersectRes:
        infoList = str(i).strip("\n").split("\t")
        transBed12 = Bed12("\t".join(infoList[:12]))
        transExonStarts, transExonEnds = transBed12.exonStarts, transBed12.exonEnds
        transExons = [(transExonStarts[x], transExonEnds[x]) for x in range(len(transExonStarts))]
        readExonStarts, readExonEnds = bedObj.reads[infoList[-3]].exonStarts, bedObj.reads[infoList[-3]].exonEnds
        readExons = [(readExonStarts[x], readExonEnds[x]) for x in range(len(readExonStarts))]
        overlapLen = getBlockLength(getOverlapOfTuple(readExons, transExons))
        coverageOnRead = overlapLen / float(getBlockLength((readExons)))
        coverageOnTrans = overlapLen / float(getBlockLength(transExons))
        if coverageOnRead >= minCoverage or coverageOnTrans >= minCoverage:
            if infoList[-3] not in singleExonStrandDict:
                singleExonStrandDict[infoList[-3]] = [transBed12.strand]
            else:
                singleExonStrandDict[infoList[-3]].append(transBed12.strand)
        else:
            singleExonStrandDict[infoList[-3]] = [bedObj.reads[infoList[-3]].strand]

    # multi Exon process
    multiExonReadDinucleotideBedObj = pybedtools.BedTool("\n".join(multiExonReadDinucleotideBedList), from_string=True)
    multiExonReadDinucleotideBedRes = multiExonReadDinucleotideBedObj.sequence(genomeFasta, name=True, tab=True)
    multiExonReadStrandDict = {}
    for i in str(open(multiExonReadDinucleotideBedRes.seqfn).read()).split("\n")[:-1]:
        infoList = str(i).strip("\n").split("\t")
        rName, dinucleotideType = infoList[0].split(":")[0], infoList[0].split(":")[2]
        if rName not in multiExonReadStrandDict:
            multiExonReadStrandDict[rName] = {"fwdCanonicalCounter": 0, "revCanonicalCounter": 0}
        if dinucleotideType == "left":
            if re.match("(G[TC])|(AT)", infoList[1]):
                multiExonReadStrandDict[rName]["fwdCanonicalCounter"] += 1
            elif re.match("[CG]T", infoList[1]):
                multiExonReadStrandDict[rName]["revCanonicalCounter"] += 1
        if dinucleotideType == "right":
            if re.match("A[GC]", infoList[1]):
                multiExonReadStrandDict[rName]["fwdCanonicalCounter"] += 1
            elif re.match("([AG]C)|(AT)", infoList[1]):
                multiExonReadStrandDict[rName]["revCanonicalCounter"] += 1

    # print
    multiExonAmbiStrandBedList = []
    for i in bedObj.reads:
        if i in singleExonStrandDict and len(set(singleExonStrandDict[i])) == 1:
            if bedObj.reads[i].strand == singleExonStrandDict[i][0]:
                print >>strandConfirmedOut, bedObj.reads[i]
            bedObj.reads[i].strand = singleExonStrandDict[i][0]
            print >>strandAdjustOut, bedObj.reads[i]
        elif i in multiExonReadStrandDict:
            if multiExonReadStrandDict[i]["fwdCanonicalCounter"] - multiExonReadStrandDict[i][
                "revCanonicalCounter"] > juncDiffScore:
                if bedObj.reads[i].strand == "+":
                    print >>strandConfirmedOut, bedObj.reads[i]
                bedObj.reads[i].strand = "+"
                print >> strandAdjustOut, bedObj.reads[i]
            elif multiExonReadStrandDict[i]["revCanonicalCounter"] - multiExonReadStrandDict[i][
                "fwdCanonicalCounter"] > juncDiffScore:
                if bedObj.reads[i].strand == "-":
                    print >>strandConfirmedOut, bedObj.reads[i]
                bedObj.reads[i].strand = "-"
                print >> strandAdjustOut, bedObj.reads[i]
            else:
                chrom, start, end, strand = bedObj.reads[i].chrom, bedObj.reads[i].chromStart, bedObj.reads[i].chromEnd, \
                                            bedObj.reads[i].strand
                multiExonAmbiStrandBedList.append(" ".join(map(str, [chrom, start, end, i, ".", strand])))
        else:
            print >>strandAdjustOut, bedObj.reads[i]
            print >>strandConfirmedOut, bedObj.reads[i]
    multiExonAmbiStrandBedObj = pybedtools.bedtool.BedTool("\n".join(multiExonAmbiStrandBedList), from_string=True)
    multiExonAmbiStrandIntersectRes = gpeBedObj.intersect(multiExonAmbiStrandBedObj, wa=True, wb=True)
    multiExonAmbiStrandDict = {}
    for i in multiExonAmbiStrandIntersectRes:
        infoList = str(i).strip("\n").split("\t")
        if infoList[-3] not in multiExonAmbiStrandDict:
            multiExonAmbiStrandDict[infoList[-3]] = [Bed12(str(i)).strand]
        else:
            multiExonAmbiStrandDict[infoList[-3]].append(Bed12(str(i)).strand)
    for i in multiExonAmbiStrandBedList:
        infoList = i.split(" ")
        if infoList[3] in multiExonAmbiStrandDict and set(multiExonAmbiStrandDict[infoList[3]]) == 1:
            if bedObj.reads[infoList[3]].strand == multiExonAmbiStrandDict[infoList[3]][0]:
                print >>strandConfirmedOut, bedObj.reads[infoList[3]]
            bedObj.reads[infoList[3]].strand = multiExonAmbiStrandDict[infoList[3]][0]
            print >>strandAdjustOut, bedObj.reads[infoList[3]]
        else:
            print >>strandAdjustOut, bedObj.reads[infoList[3]]
            print >>strandConfirmedOut, bedObj.reads[infoList[3]]
    strandAdjustOut.close()
    strandConfirmedOut.close()

def buildSpliceGrapherModel(refParams=None, paramObj=None):
    print str(datetime.datetime.now()) + " Start SpliceGrapher build classification modeling for group {}...".format(paramObj.group)
    prevDir = os.getcwd()
    resolveDir(os.path.join(refParams.out_dir, "SG_models", paramObj.projectName))
    cmd = "generate_splice_site_data.py -f {} -m {} -r spliceSiteFreq.txt 1>/dev/null".format(refParams.ref_genome,
                                                                                              refParams.ref_gtf)
    subprocess.call(cmd, shell=True)
    juncFreqDict = {"donor": [], "acceptor": []}
    with open("spliceSiteFreq.txt") as f:
        lineList = [i.strip() for i in f.readlines()]
        juncIndex = lineList.index("Breakdown of splice junctions:")
        sumFreq = 0
        for i in range(juncIndex + 1, len(lineList)):
            junc, info = re.split(":", lineList[i].strip())
            donor, acceptor = junc.strip().split("-")
            freq = float(re.search("\((.*)%\)", info).groups()[0])
            sumFreq += freq
            juncFreqDict["donor"].append(donor)
            juncFreqDict["acceptor"].append(acceptor)
            if sumFreq >= 98: break

    donorClass = ",".join(set(juncFreqDict["donor"]))
    acceptorClass = ",".join(set(juncFreqDict["acceptor"]))
    cmd = "build_classifiers.py -m {} -f {} -n 2000 -d {} -a {}".format(refParams.ref_gtf, refParams.ref_genome,
                                                                        donorClass, acceptorClass)
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " SpliceGrapher build classification modeling for group {} done!".format(paramObj.group)

# def collapse_bak(paramObj=None, refParams=None, useSqanti=False):
#     print str(datetime.datetime.now()) + " Collapse with cDNA_cupcake for project {} group {}...".format(
#         paramObj.projectName, paramObj.group)
#     prevDir = os.getcwd()
#     resolveDir("collapse")
#     makeLink("../filterByRefInfo/processed.bed12+", "processed.bed12+")
#     makeLink("../../aln.sorted.sam", "aln.sorted.sam")
#     if useSqanti:
#         makeLink("../sqanti/sqanti_artifact_removed.fa", "sqanti_artifact_removed.fa")
#         faBeforeCollapse = "sqanti_artifact_removed.fa"
#         filteredIds = getFxSequenceId("sqanti_artifact_removed.fa", isFa=True)
#         getSubSamByName("aln.sorted.sam", nameList=filteredIds, isBam=False, nameListIsFile=False, outPrefix="artifact_removed",
#                         sort=True, threads=refParams.threads)
#         filteredSortedSam = "artifact_removed.sorted.sam"
#     else:
#         makeLink("../filterByRefInfo/processed.fa", "processed.fa")
#         faBeforeCollapse = "processed.fa"
#         filteredIds = getFxSequenceId("processed.fa", isFa=True)
#         getSubSamByName("aln.sorted.sam", nameList=filteredIds, isBam=False, nameListIsFile=False, outPrefix="processed",
#                         sort=True, threads=refParams.threads)
#         filteredSortedSam = "processed.sorted.sam"
#     cmd = "collapse_isoforms_by_sam.py --input {} -s {} -i 0.9 -c 0.9 --dun-merge-5-shorter -o tofu 1>/dev/null 2>&1".format(filteredSortedSam, faBeforeCollapse)
#     subprocess.call(cmd, shell=True)
#     cmd = "gtfToGenePred tofu.collapsed.gff tofu.collapsed.gpe -genePredExt"
#     subprocess.call(cmd, shell=True)
#     cmd = "gpe2bed.pl tofu.collapsed.gpe -g > tofu.collapsed.bed12+"
#     subprocess.call(cmd, shell=True)
#     readsAssign(refParams.ref_bed, "tofu.collapsed.bed12+", readsColNum=13, outPrefix="tofu.collapsed.assigned",
#                 group=False)
#     processedIds = getFxSequenceId(fxFile="processed.fa", isFa=True)
#     filter(originFile="tofu.ignored_ids.txt", targetFile=processedIds, mode="e",
#            outFile="processed.ignored_id_removed.read.lst")
#     filter(originFile="processed.ignored_id_removed.read.lst", targetFile="processed.bed12+", targetField=4, mode="i",
#            outFile="processed.ignore_id_removed.bed12+")
#     readsAssign(refParams.ref_bed, "processed.ignore_id_removed.bed12+", readsColNum=14, outPrefix="reads.assigned",
#                 group=True)
#     cmd = "cut -f 1-12,15 reads.assigned.unambi.bed12+| bed2gpe.pl -b 12 -g 13 - | genePredToGtf file stdin reads.unambi.gtf -source=iFLAS"
#     subprocess.call(cmd, shell=True)
#     os.chdir(prevDir)
#     print str(datetime.datetime.now()) + " Collapse with cDNA_cupcake for project {} group {} done!".format(
#         paramObj.projectName, paramObj.group)

def collapse_no_sqanti(tgsSample=None, refParams=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Collapse with cDNA_cupcake for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    resolveDir("collapse")

    logDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "log")
    makeLink("../filterByRefInfo/processed.fa", "processed.fa")
    makeLink("../filterByRefInfo/processed.bed12+", "processed.bed12+")
    # makeLink("../../rawCorrection/flnc.mm2.sorted.fa", "flnc.mm2.sorted.fa")
    makeLink("../../rawCorrection/flnc.mm2.sorted.bam", "flnc.mm2.sorted.bam")
    processedIds = getFxSequenceId("processed.fa", isFa=True)
    getSubSamByName("flnc.mm2.sorted.bam", nameList=processedIds, isBam=False, nameListIsFile=False, outPrefix="processed", sort=True, threads=tgsSample.threads)
    if tgsSample.tgsPlat == "pacbio":
        cmd = "collapse_isoforms_by_sam.py --input processed.fa -s processed.sorted.sam --max_3_diff 75 --flnc_coverage 2 -c 0.9 -i 0.9 --dun-merge-5-shorter -o tofu 1>{}/tofu.collapse.log 2>&1".format(logDir)
    else:
        cmd = "collapse_isoforms_by_sam.py --input processed.fa -s processed.sorted.sam --max_3_diff 75 --flnc_coverage 5 -c 0.9 -i 0.9 --dun-merge-5-shorter -o tofu 1>{}/tofu.collapse.log 2>&1".format(logDir)
    subprocess.call(cmd, shell=True)
    cmd = "gtfToGenePred tofu.collapsed.good.gff tofu.collapsed.gpe -genePredExt"
    subprocess.call(cmd, shell=True)
    cmd = "gpe2bed.pl tofu.collapsed.gpe -g > tofu.collapsed.bed12+"
    subprocess.call(cmd, shell=True)
    readsAssign(refParams.ref_bed, "tofu.collapsed.bed12+", readsColNum=13, outPrefix="tofu.collapsed.assigned", group=True)
    # with open("tofu.collapsed.assigned.bed12+") as f:
    #     outTmp = open("tofu.collapsed.assigned.uniq.bed12+", "w")
    #     for line in f.readlines():
    #         infoList = line.strip("\n").split("\t")
    #         if infoList[13] == "IG": continue
    #         # if len(set(infoList[15].split(","))) != 1: continue
    #         print >>outTmp, "\t".join(infoList[0:13]) + "\t" + infoList[15].split(",")[0]
    #     outTmp.close()
    # filter(originFile="tofu.ignored_ids.txt", targetFile=processedIds, mode="e", outFile="processed.ignored_id_removed.read.lst")
    # cmd = "seqkit grep processed.fa -f processed.ignored_id_removed.read.lst -w 0 > processed.ignore_id_removed.fa"
    cmd = '''seqkit grep processed.fa -f <(cut -f 2 tofu.collapsed.group.txt | tr ',' '\n') -w 0 > processed.ignore_id_removed.fa'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    # filter(originFile="processed.ignored_id_removed.read.lst", targetFile="processed.bed12+", targetField=4, mode="i", outFile="processed.ignore_id_removed.bed12+")
    cmd = '''filter.pl -o <(cut -f 2 tofu.collapsed.group.txt | tr ',' '\n') processed.bed12+ -2 4 -m i > processed.ignore_id_removed.bed12+'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    readsAssign(refParams.ref_bed, "processed.ignore_id_removed.bed12+", readsColNum=14, outPrefix="reads.assigned", group=True)
    # with open("reads.assigned.bed12+") as f:
    #     outTmp = open("reads.assigned.uniq.bed12+", "w")
    #     for line in f.readlines():
    #         infoList = line.strip("\n").split("\t")
    #         if infoList[14] == "IG": continue
    #         # if len(set(infoList[16].split(","))) != 1: continue
    #         print >>outTmp, "\t".join(infoList[0:14]) + "\t" + infoList[16].split(",")[0]
    #     outTmp.close()

    cmd = "cut -f 1-12,15 reads.assigned.unambi.bed12+ | bed2gpe.pl -b 12 -g 13 - | genePredToGtf file stdin reads.unambi.gtf -source=iFLAS"
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Collapse with cDNA_cupcake for project {} entry {} done!".format(
        tgsSample.projectName, tgsSample.sampleName)


def collapse(tgsPlat="pacbio", tgsSample=None, refParams=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Collapse with cDNA_cupcake for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    resolveDir("collapse")

    logDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "log")
    makeLink("../sqanti/sqanti_artifact_removed.fa", "sqanti_artifact_removed.fa")
    makeLink("../filterByRefInfo/processed.bed12+", "processed.bed12+")
    # makeLink("../../rawCorrection/flnc.mm2.sorted.fa", "flnc.mm2.sorted.fa")
    if tgsPlat == "pacbio":
        cmd = "minimap2 -ax splice:hq -uf --secondary=no -t {} {} sqanti_artifact_removed.fa 2>{}/pacbio.sqanti.mm2.log | samtools view -bS | samtools sort | samtools view -h -O sam > artifact_removed.sorted.sam".format(tgsSample.threads, refParams.ref_mm2_index, logDir)
    else:
        cmd = "minimap2 -ax splice -uf --secondary=no -t {} -k14 {} sqanti_artifact_removed.fa 2>{}/ont.sqanti.mm2.log | samtools view -bS | samtools sort | samtools view -h -O sam > artifact_removed.sorted.sam".format(tgsSample.threads, refParams.ref_mm2_index, logDir)
    subprocess.call(cmd, shell=True)
    cmd = "collapse_isoforms_by_sam.py --input sqanti_artifact_removed.fa -s artifact_removed.sorted.sam --max_3_diff 75 --flnc_coverage 3 --dun-merge-5-shorter -o tofu 1>{}/tofu.collapse.log 2>&1".format(logDir)
    subprocess.call(cmd, shell=True)
    cmd = "gtfToGenePred tofu.collapsed.good.gff tofu.collapsed.gpe -genePredExt"
    subprocess.call(cmd, shell=True)
    cmd = "gpe2bed.pl tofu.collapsed.gpe -g > tofu.collapsed.bed12+"
    subprocess.call(cmd, shell=True)
    readsAssign(refParams.ref_bed, "tofu.collapsed.bed12+", readsColNum=13, outPrefix="tofu.collapsed.assigned", group=True)
    with open("tofu.collapsed.assigned.bed12+") as f:
        outTmp = open("tofu.collapsed.assigned.uniq.bed12+", "w")
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            if infoList[13] == "IG": continue
            if len(set(infoList[15].split(","))) != 1: continue
            print >>outTmp, "\t".join(infoList[0:13]) + "\t" + infoList[15].split(",")[0]
        outTmp.close()

    cmd = "seqkit seq -i -n sqanti_artifact_removed.fa > artifact_removed.read.lst"
    subprocess.call(cmd, shell=True)
    cmd = '''filter.pl -o <(cut -f 2 tofu.collapsed.group.txt | tr ',' '\n' ) artifact_removed.read.lst -m i |
             tee artifact_removed.ignored_id_removed.read.lst |
             seqkit grep sqanti_artifact_removed.fa -f - -w 0 > sqanti_artifact_removed.ignore_id_removed.fa
    '''
    subprocess.call(cmd, shell=True)
    cmd = "filter.pl -o artifact_removed.ignored_id_removed.read.lst processed.bed12+ -2 4 -m i > processed.ignore_id_removed.bed12+"
    subprocess.call(cmd, shell=True)
    readsAssign(refParams.ref_bed, "processed.ignore_id_removed.bed12+", readsColNum=14, outPrefix="reads.assigned",
                group=True)
    with open("reads.assigned.bed12+") as f:
        outTmp = open("reads.assigned.uniq.bed12+", "w")
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            if infoList[14] == "IG": continue
            if len(set(infoList[16].split(","))) != 1: continue
            print >>outTmp, "\t".join(infoList[0:14]) + "\t" + infoList[16].split(",")[0]
        outTmp.close()
    cmd = "cut -f 1-12,15 reads.assigned.uniq.bed12+| bed2gpe.pl -b 12 -g 13 - | genePredToGtf file stdin reads.unambi.gtf -source=iFLAS"
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Collapse with cDNA_cupcake for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)

def filterByReferenceInfo(refParams=None, tgsSample=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Filter by reference information for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    resolveDir("filterByRefInfo")
    mappedBam = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "rawCorrection",  "flnc.mm2.sorted.bam")
    mappedFa = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "rawCorrection",  "flnc.mm2.sorted.fa")
    cmd = r'''samAddTag.pl --checkHardClip --coverage --identity {} 2>lengthInconsistent.sam | samtools sort -m 4G - >mapped.addCVandID.bam'''.format(
        mappedBam)
    subprocess.call(cmd, shell=True)
    cmd = "sam2bed.pl -t CV,ID mapped.addCVandID.bam >mapped.addCVandID.bed12+"
    subprocess.call(cmd, shell=True)
    cmd = "readsFilter.pl -c 0.4 -r 0.8 mapped.addCVandID.bed12+ 2>discarded.bed12+ >uniq.bed12+"
    subprocess.call(cmd, shell=True)
    cmd = r'''(samtools view -H mapped.addCVandID.bam; samtools view mapped.addCVandID.bam | filter.pl -o <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2+1,$4,$5}' uniq.bed12+) -1 1,2,3,4 -2 3,4,1,5 -m i) > uniq.sam'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    pysam.sort("-o", "uniq.sorted.bam", "-@", str(tgsSample.threads), "uniq.sam")

    strandAdjust(refParams.ref_genome, refParams.ref_gpe, "uniq.bed12+", 2, 0.8, strandAdjust="strandAdjusted.bed12+",
                 strandConfirmed="strandConfirm.bed12+")

    if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
        if tgsSample.ngsJunctions == None:
            tgsSample.ngsJunctions = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "RNA-seq", "reassembly", "junctions.bed12")
        juncScoringParams = "-r {} strandConfirm.bed12+ -j {}".format(refParams.ref_gpe, tgsSample.ngsJunctions)
    else:
        juncScoringParams = "-r {} strandConfirm.bed12+".format(refParams.ref_gpe)
    cmd = "juncConsensus.pl -s <(juncScoring.pl {}) -l 10 strandConfirm.bed12+ >processed.bed12+".format(juncScoringParams)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    # cmd = r'''(samtools view -F 0x10 uniq.sorted.bam |
    #            perl -ne 'print if (split "\t")[5] =~ /(\d+)S$/ && $1 >30';samtools view -f 0x10 uniq.sorted.bam |
    #            perl -ne 'print if (split "\t")[5] =~ /^(\d+)S/ && $1 >30') |
    #            filter.pl -o /dev/stdin -2 4 consensus.bed12+ >processed.bed12+'''
    # subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "seqkit grep {} -f <(cut -f4 processed.bed12+) -w 0 >processed.fa".format(mappedFa)
    subprocess.call(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE)

    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Filter by reference information for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)


# def filterByJunction(refParams=None, paramObj=None):
#     print str(datetime.datetime.now()) + " Filter by junction..."
#     prevDir = os.getcwd()
#     resolveDir("filterByJunc")
#     refMD5File = os.path.join(os.path.dirname(refParams.ref_gtf), "refMD5.txt")
#     refGtfMD5 = hashlib.md5(refParams.ref_gtf).hexdigest()
#     if not os.path.exists(refMD5File):
#         md5out = open(refMD5File, "w")
#         SGModelPath = os.path.join(refParams.out_dir, "SG_models", paramObj.projectName)
#         print >>md5out, "\t".join([refParams.ref_gtf, refGtfMD5, SGModelPath])
#         md5out.close()
#         buildSpliceGrapherModel(refParams=refParams, paramObj=paramObj)
#     else:
#         with open(refMD5File) as refMD5:
#             for i in refMD5.readlines():
#                 if i.split("\t")[1] == refGtfMD5:
#                     SGModelPath = i.strip().split("\t")[2]
#                     break
#             else:
#                 md5out = open(refMD5File, "wa")
#                 SGModelPath = os.path.join(refParams.out_dir, "SG_models", paramObj.projectName)
#                 print >> md5out, "\t".join([refParams.ref_gtf, refGtfMD5, SGModelPath])
#                 md5out.close()
#                 buildSpliceGrapherModel(refParams=refParams, paramObj=paramObj)
#     mappedSam = os.path.join(refParams.out_dir, "aln.sam")
#     cmd = r'''samAddTag.pl --checkHardClip --coverage --identity --unmapped unmapped.sam {} 2>lengthInconsistent.sam | samtools view -buS - | samtools sort -m 4G - >mapped.sorted.bam'''.format(mappedSam)
#     subprocess.call(cmd, shell=True)
#     cmd = "sam2bed.pl -t CV,ID mapped.sorted.bam >mapped.bed12+"
#     subprocess.call(cmd, shell=True)
#     cmd = "readsFilter.pl -r 0.8 mapped.bed12+ 2>discarded.bed12+ >uniq.bed12+"
#     subprocess.call(cmd, shell=True)
#     cmd = r'''(samtools view -H mapped.sorted.bam; samtools view mapped.sorted.bam | filter.pl -o <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2+1,$4,$5}' uniq.bed12+) -1 1,2,3,4 -2 3,4,1,5 -m i) > uniq.sam'''
#     subprocess.call(cmd, shell=True, executable="/bin/bash")
#     print str(datetime.datetime.now()) + " Start filtering sam file with model built by SpliceGrapher..."
#     cmd = "sam_filter.py uniq.sam {}/classifiers.zip -f {} -o filtered.uniq.sam".format(SGModelPath, refParams.ref_genome)
#     subprocess.call(cmd, shell=True)
#     cmd = r'''(grep '^@' filtered.uniq.sam; grep -v '^@' filtered.uniq.sam| awk 'OFS="\t"{gsub(/^\**$/, "*", $11);print $0}') > filtered.uniq.modify.sam'''
#     subprocess.call(cmd, shell=True, executable="/bin/bash")
#     cmd = "rm *_acc.cfg *_acc.fa *_acc.svm *_don.cfg *_don.fa *_don.svm"
#     subprocess.call(cmd, shell=True)
#     # pysam.sort("-m", "4G", "-@", "5", "-o", "uniq1.sorted.bam", "filtered.uniq.sorted.sam")
#     cmd = "samtools view -bS filtered.uniq.modify.sam | samtools sort -m 4G >filtered.uniq.sorted.bam"
#     subprocess.call(cmd, shell=True)
#     print str(datetime.datetime.now()) + " Filtering sam file with model built by SpliceGrapher done!"
#
#     cmd = "strandAdjust.pl -f {} -g {} -c 0.8 -s 2 uniq.bed12+ >strandAdjusted.bed12+".format(refParams.ref_genome,
#                                                                                               refParams.ref_gpe)
#     subprocess.call(cmd, shell=True)
#
#     # confirm strand
#     uniqBed = BedFile("uniq.bed12+", type="bed12+")
#     strandAdjustedBed = BedFile("strandAdjusted.bed12+", type="bed12+")
#     strandConfirm = open("strandConfirm.bed12+", "w")
#
#     for read in uniqBed.reads:
#         if uniqBed.reads[read].strand != strandAdjustedBed.reads[read].strand: continue
#         print >> strandConfirm, uniqBed.reads[read]
#     strandConfirm.close()
#     # cmd = "paste <(cut -f4,6 uniq.bed12+) <(cut -f6 strandAdjusted.bed12+) | awk '$2!=$3{print $1}' >strandConflict.readName"
#     # cmd = "filter.pl -o strandConflict.readName -2 4 uniq.bed12+ >strandConfirm.bed12+"
#
#     # if paramObj.ngs_cfg_file:
#     #     cmd = "removeMisAlignExon.pl -g ../../raw.gpe -j $junction strandConfirm.bed12+ | fillMissExon.pl -g ../../raw.gpe -b -j $junction >exonRealign.bed12+"
#     #     cmd = '''juncScoringParams="{} -j $junction"'''.format(juncScoringParams)
#     # else:
#     if paramObj.ngsProject != None:
#         cmd = "removeMisAlignExon.pl -g {} -j {} strandConfirm.bed12+ | fillMissExon.pl -g {} -b -j /dev/null >exonRealign.bed12+".format(
#             refParams.ref_gpe, paramObj.ngsJunctions, refParams.ref_gpe)
#         subprocess.call(cmd, shell=True, executable="/bin/bash")
#         juncScoringParams = "-r {} exonRealign.bed12+ -j {}".format(refParams.ref_gpe, paramObj.ngsJunctions)
#     else:
#         cmd = "removeMisAlignExon.pl -g {} -j /dev/null strandConfirm.bed12+ | fillMissExon.pl -g {} -b -j /dev/null >exonRealign.bed12+".format(
#             refParams.ref_gpe, refParams.ref_gpe)
#         subprocess.call(cmd, shell=True, executable="/bin/bash")
#         juncScoringParams = "-r {} exonRealign.bed12+".format(refParams.ref_gpe)
#     cmd = "juncConsensus.pl -s <(juncScoring.pl {}) exonRealign.bed12+ >consensus.bed12+".format(juncScoringParams)
#     subprocess.call(cmd, shell=True, executable="/bin/bash")
#     # cmd = "bed2gpe.pl consensus.bed12+ | gpeFeature.pl --exon | filter.pl -o /dev/stdin -1 4 -2 4 consensus.bed12+ >noGap.bed12+"
#     # subprocess.call(cmd, shell=True)
#     # cmd = '''dnaOrInternalPrimingContFilter.pl -b {} -t <(sed 's/^@//' ../paTailLength.tsv) -r dnaOrInternalPrimingContFilter.log consensus.bed12+ >deCont.bed12+ 2>dnaCont.bed12+'''.format(refParams.ref2bit)
#     # subprocess.call(cmd, shell=True, executable="/bin/bash")
#     cmd = r'''(samtools view -F 0x10 filtered.uniq.sorted.bam | perl -ne 'print if (split "\t")[5] =~ /(\d+)S$/ && $1 >30';samtools view -f 0x10 filtered.uniq.sorted.bam | perl -ne 'print if (split "\t")[5] =~ /^(\d+)S/ && $1 >30') | filter.pl -o /dev/stdin -2 4 consensus.bed12+ >processed.bed12+'''
#     subprocess.call(cmd, shell=True, executable="/bin/bash")
#     cmd = "samtools fasta {} 2>/dev/null | seqkit grep - -f <(cut -f4 processed.bed12+) -w 0 >processed.fa".format(mappedSam)
#     subprocess.call(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE)
#
#     os.chdir(prevDir)
#     print str(datetime.datetime.now()) + " Filter by junction done!"

def sqantiRemoveArtifacts(refParams=None, tgsSample=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Remove artifact isoforms for project {} group {} by using SQANTI...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    resolveDir("sqanti")
    logDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "log")
    makeLink("../filterByRefInfo/processed.bed12+", "processed.bed12+")
    makeLink("../filterByRefInfo/processed.fa", "processed.fa")
    cmd = "bedToGenePred processed.bed12+ processed.gp"
    subprocess.call(cmd, shell=True)
    cmd = "genePredToGtf file processed.gp processed.gtf -source=iFLAS"
    subprocess.call(cmd, shell=True)
    cmd = "sqanti_qc.py processed.gtf {} {} -g -d ./ -o sqanti 1>{}/sqanti_qc.log 2>&1".format(refParams.ref_gtf, refParams.ref_genome, logDir)
    subprocess.call(cmd, shell=True)
    cmd = "sqanti_filter.py sqanti_classification.txt 1>{}/sqanti_filter.log 2>&1".format(logDir)
    subprocess.call(cmd, shell=True)
    cmd = "seqkit grep processed.fa -w 0 -f Filter_out/sqanti_classification.txt_curatedTranscriptome.txt >sqanti_artifact_removed.fa"
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Remove artifact isoforms for project {} group {} by using SQANTI done!".format(tgsSample.projectName, tgsSample.sampleName)
