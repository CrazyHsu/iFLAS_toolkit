#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: refine.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-05-25 15:14:26
Last modified: 2021-05-25 15:14:26
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
            if multiExonReadStrandDict[i]["fwdCanonicalCounter"] - multiExonReadStrandDict[i]["revCanonicalCounter"] > juncDiffScore:
                if bedObj.reads[i].strand == "+":
                    print >>strandConfirmedOut, bedObj.reads[i]
                bedObj.reads[i].strand = "+"
                print >> strandAdjustOut, bedObj.reads[i]
            elif multiExonReadStrandDict[i]["revCanonicalCounter"] - multiExonReadStrandDict[i]["fwdCanonicalCounter"] > juncDiffScore:
                if bedObj.reads[i].strand == "-":
                    print >>strandConfirmedOut, bedObj.reads[i]
                bedObj.reads[i].strand = "-"
                print >> strandAdjustOut, bedObj.reads[i]
            else:
                chrom, start, end, strand = bedObj.reads[i].chrom, bedObj.reads[i].chromStart, bedObj.reads[i].chromEnd, bedObj.reads[i].strand
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

def identifyNovelIsoformsByJunctions(gpeFile, bedFile, anno="annoIsoform.bed", novel="novelIsoform.bed"):
    annoOut = open(anno, "w")
    novelOut = open(novel, "w")
    gpeObj = GenePredObj(gpeFile, bincolumn=False)
    bedObj = BedFile(bedFile, type="bed12+")
    for i in bedObj.reads:
        bedGene = bedObj.reads[i].otherList[0]
        readsIntrons = bedObj.reads[i].introns
        if len(readsIntrons) < 1: continue
        if bedGene in gpeObj.geneName2gpeObj:
            for trans in gpeObj.geneName2gpeObj[bedGene]:
                transIntrons = trans.introns
                if not len(set(readsIntrons) - set(transIntrons)):
                    print >> annoOut, bedObj.reads[i]
                    break
            else:
                print >> novelOut, bedObj.reads[i]
        else:
            print >> novelOut, bedObj.reads[i]
    annoOut.close()
    novelOut.close()

def refineJunc(dataObj=None, refParams=None, dirSpec=None, refine=True):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Refine the collapsed isoforms for project {} sample {}...".format(projectName, sampleName)
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    refineDir = os.path.join(baseDir, "refine")
    resolveDir(refineDir)
    processedFa = os.path.join(baseDir, "mapping", "flnc.processed.fa")
    processedBed = os.path.join(baseDir, "mapping", "flnc.addCVandID.bed12+")

    collapsedGff = os.path.join(baseDir, "collapse", "tofu.collapsed.good.gff")
    collapsedGroup = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    cmd = "gtfToGenePred {} tofu.collapsed.gpe -genePredExt".format(collapsedGff)
    subprocess.call(cmd, shell=True)
    cmd = "gpe2bed.pl tofu.collapsed.gpe -g > tofu.collapsed.bed12+"
    subprocess.call(cmd, shell=True)

    if refine:
        strandAdjust(refParams.ref_genome, refParams.ref_gpe, "tofu.collapsed.bed12+", 0.8, 2,
                     strandAdjust="tofu.strandAdjusted.bed12+", strandConfirmed="tofu.strandConfirm.bed12+")
        if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
            if dataObj.ngs_junctions == None:
                dataObj.ngs_junctions = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "junctions.bed")
            juncScoringParams = "-r {} tofu.strandConfirm.bed12+ -j {}".format(refParams.ref_gpe, dataObj.ngs_junctions)
        else:
            juncScoringParams = "-r {} tofu.strandConfirm.bed12+".format(refParams.ref_gpe)
        cmd = "juncConsensus.pl -s <(juncScoring.pl {}) -l 5 tofu.strandConfirm.bed12+ > tofu.juncAdjusted.bed12+".format(juncScoringParams)
        subprocess.call(cmd, shell=True, executable="/bin/bash")

        readsAssign(refParams.ref_bed, "tofu.juncAdjusted.bed12+", readsColNum=13, outPrefix="tofu.collapsed.assigned", group=True)
        # cmd = "cut -f1-12,14 tofu.collapsed.assigned.unambi.bed12+ > isoformGrouped.bed12+"
        # subprocess.call(cmd, shell=True)
    else:
        readsAssign(refParams.ref_bed, "tofu.collapsed.bed12+", readsColNum=13, outPrefix="tofu.collapsed.assigned", group=True)
    cmd = "cut -f1-12,14 tofu.collapsed.assigned.unambi.bed12+ > isoformGrouped.bed12+"
    subprocess.call(cmd, shell=True)
    identifyNovelIsoformsByJunctions(refParams.ref_gpe, "isoformGrouped.bed12+", anno="isoformGrouped.anno.bed12+",
                                     novel="isoformGrouped.novel.bed12+")

    ######## for reads
    cmd = '''seqkit grep {} -f <(cut -f 2 {} | tr ',' '\n') -w 0 > flnc.processed.ignore_id_removed.fa'''.format(processedFa, collapsedGroup)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''filter.pl -o <(cut -f 2 {} | tr ',' '\n') {} -2 4 -m i > flnc.processed.ignore_id_removed.bed12+'''.format(collapsedGroup, processedBed)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    readsAssign(refParams.ref_bed, "flnc.processed.ignore_id_removed.bed12+", readsColNum=14, outPrefix="reads.assigned", group=True)

    cmd = "cut -f 1-12,15 reads.assigned.unambi.bed12+ | bed2gpe.pl -b 12 -g 13 - | genePredToGtf file stdin reads.unambi.gtf -source=iFLAS"
    subprocess.call(cmd, shell=True)
    print getCurrentTime() + " Refine the collapsed isoforms for project {} sample {} done!".format(projectName, sampleName)
