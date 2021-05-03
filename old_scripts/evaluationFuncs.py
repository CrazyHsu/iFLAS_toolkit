#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: evaluationFuncs.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-08 19:09:28
Last modified: 2019-12-08 19:09:36
'''
import datetime, random, pybedtools
from Bio import SeqIO
from commonObjs import *

def checkNovelJuncAndIsoforms(gpeFile, bedFile, anno=None, novel=None):
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

def contamination(assignBed, refParams=None, contaminationOut="contamination.log"):
    readFeatureDict = {"exonic": 0, "intronic": 0, "intergenic": 0}
    readCount = 0
    with open(assignBed) as f:
        for line in f:
            readCount += 1
            infoList = line.strip("\n").split("\t")
            if infoList[14] == "E":
                readFeatureDict["exonic"] += 1
            elif infoList[14] == "I":
                readFeatureDict["intronic"] += 1
            elif infoList[14] == "IG":
                readFeatureDict["intergenic"] += 1
    refGPE = refParams.ref_gpe
    refSize = refParams.ref_size
    refGPEobj = GenePredObj(refGPE)
    exomeLength = refGPEobj.getGeneExonLength()
    intromeLength = refGPEobj.getGeneIntronLength()
    genomeLength = 0
    with open(refSize) as f:
        for line in f:
            genomeLength += int(line.strip("\n").split("\t")[1])
    intergenomeLength = genomeLength - exomeLength - intromeLength
    exomeRPKM = float(readFeatureDict["exonic"])/(exomeLength/1000)/(float(readCount)/1000000)
    intromeRPKM = float(readFeatureDict["intronic"])/(intromeLength/1000)/(float(readCount)/1000000)
    intergenomeRPKM = float(readFeatureDict["intergenic"])/(intergenomeLength/1000)/(float(readCount)/1000000)
    log = open(contaminationOut, "w")
    print >>log, "Exome\t{}".format(exomeRPKM)
    print >>log, "Introme\t{}".format(intromeRPKM)
    print >>log, "Intergenome\t{}".format(intergenomeRPKM)
    log.close()

def gcAcrossRead(fastaFile, outFile, interval=20):
    import math
    from Bio.SeqUtils import GC
    out = open(outFile, "w")
    for seq in SeqIO.parse(fastaFile, "fasta"):
        if len(seq) < interval:
            chunkSize = 1
        else:
            chunkSize = int(math.ceil(len(seq)/interval))
        gcList = [seq.name]
        for i in xrange(0, len(seq), chunkSize):
            gcList.append(GC(seq[i:i+chunkSize].seq))
        print >>out, "\t".join(map(str, gcList))
    out.close()

def preCorrEval(tgsSample=None, dirSpec=None):
    # Pre-mapping evaluation
    print str(datetime.datetime.now()) + " Pre-mapping evaluation for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    preCorrDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "preCorrEval")
    resolveDir(preCorrDir)
    if tgsSample.tgsPlat.lower() == "pacbio":
        # flncReport = os.path.join(dirSpec.out_dir, tgsSample.projectName, "retrieveData", "pacbioDemultiplexed", tgsSample.uniqName + ".flnc.report.csv")
        flncFa = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "rawCorrection", "flnc.fa")
        # cmd = "sed '1d' {} | cut -d',' -f 5 | hist.R -x='Tail Length' -y=Density -b=1 -d -x1=0 -x2=100 -p=paTailLength.pdf 2>/dev/null".format(flncReport)
        # subprocess.call(cmd, shell=True)

        cmd = "seqkit fx2tab -n -g {} | cut -f 1,2 > GC_of_raw_flnc_reads.log".format(flncFa)
        subprocess.call(cmd, shell=True)
        cmd = '''cut -f2 GC_of_raw_flnc_reads.log | distrCurve.R -d -m='GC Content of raw flnc Reads' -x='Binned GC%' -y='Fraction of Reads' -v=50 -p=GC_of_raw_flnc_reads.pdf 2>/dev/null'''
        subprocess.call(cmd, shell=True)

        gcAcrossRead(flncFa, "GC_across_raw_flnc_read.log")
        cmd = '''cut -f2- GC_across_raw_flnc_read.log | box.R -stack -nJ -ho=50 -m='GC Content across raw flnc Reads' -x=Interval -y=GC% -oS=0.5 -w=11 -p=GC_across_raw_flnc_read.pdf 2>/dev/null'''
        subprocess.call(cmd, shell=True)
    else:
        nanoporeRawSeq = os.path.join(dirSpec.out_dir, tgsSample.projectName, "retrieveData", "nanoporeRetrieve", tgsSample.sampleName, "nanoporeRawSeq.fq")
        if checkFast5Files(tgsSample.dataLocation) == "fast5":
            nanoporePaReport = os.path.join(dirSpec.out_dir, tgsSample.projectName, "retrieveData", "nanoporeRetrieve", tgsSample.sampleName, "polya_results.tsv")
            cmd = "cut -f 9 {} | sed '1d' | hist.R -x='Tail Length' -y=Density -b=1 -d -x1=0 -x2=100 -p=paTailLength.pdf 2>/dev/null".format(nanoporePaReport)
            subprocess.call(cmd, shell=True)

        cmd = "seqkit seq -i --rna2dna {} | tee nanopore.raw.fa | seqkit fx2tab -n -g | cut -f 1,2 > GC_of_raw_nanopore_reads.log".format(nanoporeRawSeq)
        subprocess.call(cmd, shell=True)
        cmd = '''cut -f2 GC_of_raw_nanopore_reads.log | distrCurve.R -d -m='GC Content of raw nanopore Reads' -x='Binned GC%' -y='Fraction of Reads' -v=50 -p=GC_of_raw_nanopore_reads.pdf 2>/dev/null'''
        subprocess.call(cmd, shell=True)

        gcAcrossRead("nanopore.raw.fa", "GC_across_raw_nanopore_read.log")
        cmd = '''cut -f2- GC_across_raw_nanopore_read.log | box.R -stack -nJ -ho=50 -m='GC Content across raw nanopore Reads' -x=Interval -y=GC% -oS=0.5 -w=11 -p=GC_across_raw_nanopore_read.pdf 2>/dev/null'''
        subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Pre-mapping evaluation for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)

def postCorrEval(refParams=None, tgsSample=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Post-mapping evaluation for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    postCorrDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "postCorrEval")
    resolveDir(postCorrDir)
    if not os.path.exists("length"):
        os.makedirs("length")
    makeLink("../rawCorrection/flnc.fa", "flnc.fa")
    makeLink("../rawCorrection/flnc.mm2.sorted.bed12", "flnc.mm2.sorted.bed12")
    makeLink("../filtration/collapse/reads.assigned.unambi.bed12+", "reads.assigned.unambi.bed12+")
    makeLink("../filtration/collapse/reads.assigned.bed12+", "reads.assigned.bed12+")
    makeLink("../filtration/collapse/tofu.collapsed.assigned.unambi.bed12+", "isoforms.assigned.unambi.bed12+")
    referenceGeneLenOut = open("length/referenceGene.lst", "w")
    referenceGeneToLenOut = open("length/referenceGeneToLen.txt", "w")
    unambiReadsLenOut = open("length/unambiReads.lst", "w")
    unambiReadsToLenOut = open("length/unambiReadsToLen.txt", "w")
    cmd = "seqkit fx2tab -l flnc.fa | cut -f 1,4 | tee length/rawFlncReadsToLen.txt | cut -f 2 > length/rawFlncReads.lst"
    subprocess.call(cmd, shell=True)
    with open(refParams.ref_gpe) as f:
        for i in f.readlines():
            print >> referenceGeneLenOut, sum(GenePredExtLine(i, bincolumn=False).blockSizes)
            print >> referenceGeneToLenOut, "\t".join(map(str, [GenePredExtLine(i, bincolumn=False).transName, sum(GenePredExtLine(i, bincolumn=False).blockSizes)]))
    with open("reads.assigned.unambi.bed12+") as f:
        for i in f.readlines():
            print >> unambiReadsLenOut, sum(Bed12Plus(i).blockSizes)
            print >> unambiReadsToLenOut, "\t".join(map(str, [Bed12Plus(i).name, sum(Bed12Plus(i).blockSizes)]))
    referenceGeneLenOut.close()
    referenceGeneToLenOut.close()
    unambiReadsLenOut.close()
    unambiReadsToLenOut.close()

    cmd = '''distrCurves.R -x1=0 -x2=10000 -d -x='Binned Length (limited in 0-10000)' -w=15 length/*.lst -b=150 -p=LengthDistribution.curve.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "boxes.R -ng -no length/*.lst -p=LengthDistribution.box.pdf 2>/dev/null"
    subprocess.call(cmd, shell=True)
    cmd = "geneCoverage.pl -g {} reads.assigned.unambi.bed12+ >coverage.log".format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True)
    cmd = '''awk '$3<20000' coverage.log | select.pl --index 3,2 | binBox.R -b=10 -w=15 -m='Coverage across Genes' -x='RefSeq Gene Length (<20000)' -y=Coverage -p=coverage.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    # Distance to TSS/TTS across Reads with Different Length
    cmd = '''distanceToEnd.pl -g {} -f length/unambiReadsToLen.txt reads.assigned.unambi.bed12+ | awk '$2!="NA"' >reads.distance2Ends.log'''.format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''select.pl --index 4,2 reads.distance2Ends.log | binBox.R -ng -no -b=5 -w=15 -m='Distance to Gene TSS' -x='Binned Read Length' -y=Distance -p=reads.distance2TSS.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''select.pl --index 4,3 reads.distance2Ends.log | binBox.R -ng -no -b=5 -w=15 -m='Distance to Gene TTS' -x='Binned Read Length' -y=Distance -p=reads.distance2TTS.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "cut -f2 reads.distance2Ends.log >reads.distanceToTSS.dis && cut -f3 reads.distance2Ends.log >reads.distanceToTTS.dis"
    subprocess.call(cmd, shell=True)
    cmd = '''distrCurve.R <reads.distanceToTSS.dis -d -x='Distance (limited in -500 to 500)' -y='Density of Reads' -m='Distance to TSS' -x1=-500 -x2=500 -b=10 -p=distance2EndOfTSS.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    cmd = '''distrCurve.R <reads.distanceToTTS.dis -d -x='Distance (limited in -500 to 500)' -y='Density of Reads' -m='Distance to TTS' -x1=-500 -x2=500 -b=10 -p=distance2EndOfTTS.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    cmd = "rm *.dis"
    subprocess.call(cmd, shell=True)

    cmd = '''cut -f1-12,14 isoforms.assigned.unambi.bed12+ > isoformGrouped.bed12+'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    checkNovelJuncAndIsoforms(refParams.ref_gpe, "isoformGrouped.bed12+", anno="isoformGrouped.anno.bed12+", novel="isoformGrouped.novel.bed12+")
    contamination("reads.assigned.bed12+", refParams=refParams, contaminationOut="contamination.log")
    cmd = '''bar.R <contamination.log -anno -c=darkgreen -f=white -m=Contamination -x='Genomic Context' -y="Abundance" -p=contamination.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)

    # Saturation of Sequencing
    saturation("flnc.mm2.sorted.bed12", refParams=refParams, tgsSample=tgsSample)
    cmd = '''lines.R *.discover -m='Saturation of Sequencing' -x='Reads Count' -y='Discovery Rate' -w=15 -lgPos=top -p=discoveryRate.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)

    # Sequencing Error Rate
    errorRate = open("readsErrorRate.log", "w")
    with open("reads.assigned.unambi.bed12+") as f:
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            print >>errorRate, "{}\t{}".format(lineInfo[3], 100-float(lineInfo[12])*float(lineInfo[13])/100)
    errorRate.close()
    cmd = '''cut -f2 readsErrorRate.log | distrCurve.R - -d -m='Rate of Mismatch & Indel for Mapping' -x=Rate -y='Fraction of reads' -p=readsErrorRate.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    # Read Count per Gene
    cmd = "geneRPKM.pl -g {} reads.assigned.unambi.bed12+ | tee RPKM.bed6+ | cut -f7 >readsCount.lst".format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True)
    cmd = '''distrCurve.R <readsCount.lst -m='Gene Count at reads Count' -x='Reads Count' -y='Gene Count' -xl=10 -yl=10 -ng -p=readsCountPerGene.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
        if tgsSample.ngsJunctions == None:
            tgsSample.ngsJunctions = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "RNA-seq", "reassembly", "junctions.bed12")
            cmd = '''awk '$10>1' isoforms.assigned.unambi.bed12+ | bed2gpe.pl | transSupportByJunction.pl -j {} >supportedByRNAseq.tsv 2>supportedByRNAseq.summary'''.format(tgsSample.ngsJunctions)
            subprocess.call(cmd, shell=True, executable="/bin/bash")
            cmd = '''awk 'BEGIN{OFS="\t"}{print $1,$2,$4/$3}' supportedByRNAseq.summary | summary2d.R -binWidthX=0.9999999 -binWidthY=0.9999999 -x='Junction Count of PacBio Reads' -y='Supported Junction Count' -fL=gray -fH=red -w=12 -p=supportedByRNAseq.pdf'''
            subprocess.call(cmd, shell=True, executable="/bin/bash")

    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Post-mapping evaluation for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)

def parallelSaturation(params):
    bedLines, count, assignBedDict, refExonBedObj = params
    random.shuffle(bedLines)
    headIreads = bedLines[0:count]
    shufReads = [b.strip().split("\t")[3] for b in headIreads]
    transReads = set()
    geneReads = set()
    for j in shufReads:
        if j in assignBedDict and assignBedDict[j][0] == "E":
            transReads.add(assignBedDict[j][1])
            geneReads.add(assignBedDict[j][2])
    shufExonBedList = []
    for l in headIreads:
        tmpRead = Bed12Plus(l.strip("\n"))
        for z in range(len(tmpRead.exonStarts)):
            shufExonBedList.append("\t".join(map(str, [tmpRead.chrom, tmpRead.exonStarts[z], tmpRead.exonEnds[z]])))
    shufExonBedList = set(shufExonBedList)
    shufExonBedObj = pybedtools.BedTool("\n".join(shufExonBedList), from_string=True)
    intersectRes = refExonBedObj.intersect(shufExonBedObj, wa=True, nonamecheck=True)
    return count, len(transReads), len(geneReads), len(set(intersectRes))

def saturation(bedFile, refParams=None, tgsSample=None):
    cmd = "readsAssigner.pl -g {} {} >flncReads.assign.multLine.bed12+".format(refParams.ref_gpe, bedFile)
    subprocess.call(cmd, shell=True)
    refExonBedList = []
    with open(refParams.ref_gpe) as f:
        for line in f.readlines():
            lineGPE = GenePredExtLine(line, bincolumn=False)
            for i in range(len(lineGPE.exonStarts)):
                refExonBedList.append("\t".join(map(str, [lineGPE.chrom, lineGPE.exonStarts[i], lineGPE.exonEnds[i]])))
    refExonBedList = set(refExonBedList)
    refExonNum = len(refExonBedList)
    refExonBedObj = pybedtools.BedTool("\n".join(refExonBedList), from_string=True)
    assignBedDict = {}
    with open("flncReads.assign.multLine.bed12+") as f:
        for line in f:
            infoList = line.strip().split("\t")
            assignBedDict[infoList[3]] = infoList[12:]
    bedLines = open(bedFile).readlines()
    totalReadN = len(bedLines)
    refGenes = [[1, line.split("\t")[11]] for line in open(refParams.ref_gpe)]
    transN = sum(i[0] for i in refGenes)
    geneN = len(set([i[1] for i in refGenes]))
    step1 = totalReadN/80
    step2 = totalReadN/40
    halfOfTotal = totalReadN/2
    transDis = open("trans.discover", "w")
    geneDis = open("gene.discover", "w")
    exonDis = open("exon.discover", "w")
    p = Pool(tgsSample.threads)
    paramsList = []
    for i in xrange(0, halfOfTotal, step1):
        paramsList.append([bedLines, i, assignBedDict, refExonBedObj])
    for i in xrange(halfOfTotal, totalReadN, step2):
        paramsList.append([bedLines, i, assignBedDict, refExonBedObj])
    result = p.map_async(parallelSaturation, paramsList)
    result.wait()
    if result.ready():
        if result.successful():
            resList = result.get()
            sorted(resList, key=lambda s: s[0])
            for j in resList:
                print >> transDis, "{}\t{}".format(j[0], j[1] / float(transN))
                print >> geneDis, "{}\t{}".format(j[0], j[2] / float(geneN))
                print >> exonDis, "{}\t{}".format(j[0], j[3] / float(refExonNum))
        else:
            print "There maybe something wrong, please check it!"
    transDis.close()
    geneDis.close()
    exonDis.close()
