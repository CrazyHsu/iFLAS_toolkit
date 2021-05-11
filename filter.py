#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: filter.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:15:41
Last modified: 2021-04-29 16:15:41
'''

from commonFuncs import *
from commonObjs import *
import pybedtools

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


def getJuncFromRegtools(dataObj=None, dirSpec=None, filterByCount=10):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Get junctions from RNA-seq with Regtools for project {} sample {}...".format(projectName, sampleName)
    rnaseqSortedBam = os.path.join(dirSpec.out_dir, projectName, sampleName, "mapping", "rna-seq", "reassembly", "tmp.bam")
    cmd = "regtools junctions extract -a 5 -m 50 -M 50000 {} > tmp.bed".format(rnaseqSortedBam)
    subprocess.call(cmd, shell=True)
    cmd = '''awk '{if($5>''' + str(filterByCount) + '''){print}}' tmp.bed > junctions.bed'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    dataObj.ngs_junctions = os.path.join(os.getcwd(), "junctions.bed")
    print getCurrentTime() + " Get junctions from RNA-seq with Regtools for project {} sample {} done!".format(projectName, sampleName)

def filterByJunc(dataObj=None, refParams=None, dirSpec=None, threads=10):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Filter reads by junctions information from RNA-seq for project {} sample {}...".format(projectName, sampleName)
    mappedSam = os.path.join(dataObj.out_dir, projectName, sampleName, "mapping", "flnc.mm2.sam")
    rawMappedSam = os.path.join(dataObj.out_dir, projectName, sampleName, "mapping", "rawFlnc.mm2.sam")
    cmd = r'''samAddTag.pl --checkHardClip --coverage --identity {} 2>lengthInconsistent.sam | samtools sort -m 4G - >mapped.addCVandID.bam'''.format(mappedSam)
    subprocess.call(cmd, shell=True)
    cmd = "sam2bed.pl -t CV,ID mapped.addCVandID.bam >mapped.addCVandID.bed12+"
    subprocess.call(cmd, shell=True)
    cmd = r'''samAddTag.pl --checkHardClip --coverage --identity {} 2>/dev/null | samtools sort -m 4G - >raw.mapped.addCVandID.bam'''.format(rawMappedSam)
    subprocess.call(cmd, shell=True)
    cmd = "sam2bed.pl -t CV,ID raw.mapped.addCVandID.bam >raw.mapped.addCVandID.bed12+"
    subprocess.call(cmd, shell=True)

    cmd = "readsFilter.pl -c 0.4 -r 0.8 mapped.addCVandID.bed12+ 2>discarded.bed12+ >uniq.bed12+"
    subprocess.call(cmd, shell=True)
    cmd = r'''(samtools view -H mapped.addCVandID.bam; samtools view mapped.addCVandID.bam | filter.pl -o <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2+1,$4,$5}' uniq.bed12+) -1 1,2,3,4 -2 3,4,1,5 -m i) > uniq.sam'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "samtools sort -@ {} uniq.sam > uniq.sorted.bam".format(threads)
    subprocess.call(cmd, shell=True)

    strandAdjust(refParams.ref_genome, refParams.ref_gpe, "uniq.bed12+", 2, 0.8, strandAdjust="strandAdjusted.bed12+",
                 strandConfirmed="strandConfirm.bed12+")

    if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
        if dataObj.ngs_junctions == None:
            dataObj.ngs_junctions = os.path.join(dirSpec.out_dir, projectName, sampleName, "RNA-seq", "reassembly", "junctions.bed")
        juncScoringParams = "-r {} strandConfirm.bed12+ -j {}".format(refParams.ref_gpe, dataObj.ngs_junctions)
    else:
        juncScoringParams = "-r {} strandConfirm.bed12+".format(refParams.ref_gpe)
    cmd = "juncConsensus.pl -s <(juncScoring.pl {}) -l 10 strandConfirm.bed12+ >processed.bed12+".format(juncScoringParams)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    flncFx = dataObj.data_processed_location
    cmd = "seqkit grep {} -f <(cut -f4 processed.bed12+) | seqkit fq2fa - -w 0 >processed.fa".format(flncFx)
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    print getCurrentTime() + " Filter reads by junctions information from RNA-seq for project {} sample {} done!".format(projectName, sampleName)

def filter(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    workDir = os.path.join(dataObj.out_dir, projectName, sampleName, "filtration")
    resolveDir(workDir)
    if dataObj.ngs_junctions == None and (dataObj.ngs_left_reads != None or dataObj.ngs_right_reads != None):
        getJuncFromRegtools(dataObj=dataObj, dirSpec=dirSpec)
    filterByJunc(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, threads=dataObj.single_run_threads)

