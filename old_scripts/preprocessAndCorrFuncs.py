#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: preprocessAndCorrFunc.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-08 17:27:39
Last modified: 2019-12-08 17:27:42
'''
import os, subprocess, datetime, shutil, sys, pysam, pybedtools, re
from commonFuncs import *

def correctWithDifferentTools(seqNeedCorrect=None, inputBam=None, paramObj=None, refParams=None, mecatParams=None,
                              useProovread=None, projectsDict=None,
                              tgsPlat="pacbio"):
    if not os.path.exists(refParams.ref_mm2_index):
        cmd = "minimap2 -d {} -t {} {}".format(refParams.ref_mm2_index, refParams.threads, refParams.ref_genome)
        subprocess.call(cmd, shell=True)
    if tgsPlat == "pacbio" and paramObj.ngsProject:
        isoseq3Corr(inputBam, paramObj=paramObj, refParams=refParams, useProovread=useProovread,
                    projectsDict=projectsDict)
    elif tgsPlat == "pacbio" and paramObj.ngsProject == None:
        isoseq3Corr(inputBam, paramObj=paramObj, refParams=refParams)
    elif tgsPlat == "nanopore" and paramObj.ngsProject:
        mecatCorr(seqNeedCorrect, paramObj=paramObj, refParams=refParams, mecatParams=mecatParams,
                  useProovread=useProovread, projectsDict=projectsDict)
    elif tgsPlat == "nanopore" and paramObj.ngsProject == None:
        mecatCorr(seqNeedCorrect, paramObj=paramObj, refParams=refParams, mecatParams=mecatParams)
    elif paramObj.ngsProject:
        print str(datetime.datetime.now()) + " No TGS reads provided, only NGS analysis will be conducted!"

    cmd = "(samtools view -H aln.sam;samtools view -f 16 -F 4079 aln.sam;samtools view -f 0 -F 4095 aln.sam) | samtools sort -@ {} --output-fmt=SAM > aln.sorted.sam".format(refParams.threads)
    subprocess.call(cmd, shell=True)
    pysam.view("-o", "aln.sorted.bam", "-@", str(refParams.threads), "--output-fmt", "BAM", "aln.sorted.sam", catch_stdout=False)
    cmd = "samtools index aln.sorted.bam"
    subprocess.call(cmd, shell=True)
    # pysam.index("aln.sorted.bam", catch_stdout=False)
    samProcess("aln.sorted.bam", isBam=True, outPrefix="aln.sorted", toFa=True)
    pysam.view("-o", "unmapped.sam", "-h", "-f", "4", "aln.sam", catch_stdout=False)
    samProcess("unmapped.sam", outPrefix="unmapped", toFa=True)
    pybedtools.BedTool("aln.sorted.bam").bam_to_bed(split=True, bed12=True).saveas("aln.sorted.bed12")

def isoseq3Corr(inputBam, tgsSample=None, refParams=None, useProovread=False, projectsDict=None):
    print str(datetime.datetime.now()) + " IsoSeq3 correction for project {} entry {}...".format(paramObj.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    inputBam = os.path.abspath(os.path.join(prevDir, inputBam))
    resolveDir("isoseq3")
    if not os.path.exists("primer"):
        os.makedirs("primer")
    logDir = os.path.join(refParams.out_dir, tgsSample.projectName, tgsSample.sampleName, "log")
    shutil.copy(paramObj.tgsPrimer, "primer/primer.fasta")
    isoseqPrimer = "primerIsoseq.fa"
    limaPrimer = "primerLima.fa"
    makePrimerForLimaAndIsoseq("primer/primer.fasta", isoseqPrimer, limaPrimer)
    cmd = "lima --isoseq --dump-clips -j {} {} primerLima.fa demux.bam".format(refParams.threads, inputBam)
    # print cmd
    subprocess.call(cmd, shell=True)

    cmd = "isoseq3 refine -j {} --require-polya demux.5p--3p.bam primerLima.fa flnc.bam".format(refParams.threads)
    # print cmd
    subprocess.call(cmd, shell=True)

    samProcess("flnc.bam", isBam=True, outPrefix="flnc", toFq=True, toFa=True)
    if useProovread:
        # hybrid strategy will use proovread to correct
        flncFq = os.path.join(os.getcwd(), "flnc.fq")
        proovreadCorr(seqNeedCorrect=flncFq, paramObj=paramObj, refParams=refParams, projectsDict=projectsDict)
        proovreadFq = os.path.join(os.getcwd(), "proovread", "proovread.fq")
        processedFlncFq = proovreadFq

        cmd = "minimap2 -ax splice:hq -uf {} -t {} --secondary=no --MD {} >aln.sam 2>{}/{}.proovreadFlnc.mm2.log".format(
            refParams.ref_mm2_index, tgsSample.threads, processedFlncFq, logDir, tgsSample.sampleName)
        subprocess.call(cmd, shell=True)
        pysam.sort("-o", "aln.raw.sorted.sam", "--output-fmt", "SAM", "-@", str(refParams.threads), "aln.sam", catch_stdout=False)
        cmd = "minimap2 -ax splice:hq -uf {} -t {} --secondary=no --MD flnc.fa >flnc.sam 2>{}/{}.rawFlnc.mm2.log".format(
            refParams.ref_mm2_index, tgsSample.threads, logDir, tgsSample.sampleName)
        subprocess.call(cmd, shell=True)
        flncSam = "flnc.sam"
    else:
        processedFlncFq = "flnc.fq"
        cmd = "minimap2 -ax splice:hq -uf {} -t {} --secondary=no --MD {} >aln.sam 2>{}/{}.rawFlnc.mm2.log".format(
            refParams.ref_mm2_index, refParams.threads, processedFlncFq, logDir, paramObj.group)
        subprocess.call(cmd, shell=True)
        pysam.sort("-o", "aln.raw.sorted.sam", "--output-fmt", "SAM", "-@", str(refParams.threads), "aln.sam", catch_stdout=False)
        flncSam = "aln.sam"
    cmd = "(samtools view -H {};samtools view -f 16 -F 4079 {};samtools view -f 0 -F 4095 {}) | samtools sort -@ {} > flnc.mm2.sorted.bam".format(
        flncSam, flncSam, flncSam, refParams.threads)
    subprocess.call(cmd, shell=True)
    pybedtools.BedTool("flnc.mm2.sorted.bam").bam_to_bed(bed12=True).saveas("flnc.mm2.sorted.bed12")

    os.chdir(prevDir)
    if os.path.islink("aln.sam") or os.path.exists("aln.sam"):
        os.remove("aln.sam")
    mappedSam = os.path.join(prevDir, "isoseq3", "aln.sam")
    makeLink(mappedSam, "aln.sam")
    print str(datetime.datetime.now()) + " IsoSeq3 correction for project {} group {} done!".format(paramObj.projectName, paramObj.group)

def makePrimerForLimaAndIsoseq(primer, isoseqPrimer, limaPrimer):
    records = list(SeqIO.parse(primer, "fasta"))
    if len(records) < 2:
        sys.exit("The primer is not correct!")
    isoseqOut = open(isoseqPrimer, "w")
    limaOut = open(limaPrimer, "w")
    for i in range(len(records)):
        if i == 0:
            print >> isoseqOut, ">F0\n{}".format(records[i].seq)
            print >> limaOut, ">5p\n{}".format(records[i].seq)
        else:
            print >> isoseqOut, ">R0\n{}".format(records[i].seq)
            print >> limaOut, ">3p\n{}".format(records[i].seq)
    isoseqOut.close()
    limaOut.close()


def isoseq3Corr_reads(inputBam, subreadsXML, paramObj=None, refParams=None, hybrid=False, projectsDict=None):
    print str(datetime.datetime.now()) + " IsoSeq3 correction..."
    prevDir = os.getcwd()
    resolveDir("isoseq3")
    if not os.path.exists("primer"):
        os.makedirs("primer")
    logDir = os.path.join(refParams.out_dir, paramObj.projectName, paramObj.group, "log")
    shutil.copy(paramObj.tgsPrimer, "primer/primer.fasta")
    isoseqPrimer = "primerIsoseq.fa"
    limaPrimer = "primerLima.fa"
    makePrimerForLimaAndIsoseq("primer/primer.fasta", isoseqPrimer, limaPrimer)
    cmd = "lima --isoseq --dump-clips -j {} ../{} primerLima.fa demux.bam".format(refParams.threads, inputBam)
    # print cmd
    subprocess.call(cmd, shell=True)

    cmd = "isoseq3 refine --require-polya demux.5p--3p.bam primerLima.fa flnc.bam"
    # print cmd
    subprocess.call(cmd, shell=True)

    cmd = '''sed '1d' flnc.report.csv | awk 'BEGIN{FS=",";OFS="\t"}{print $1,"Z:"$1}'> flnc.reads.imTag.lst'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    cmd = "samAddTag.pl --listTag=im --list=flnc.reads.lst flnc.bam | samtools view -bS > flnc.reads.imAdd.bam"
    subprocess.call(cmd, shell=True)

    cmd = "isoseq3 polish -j {} flnc.reads.imAdd.bam ../{} polished.flnc.bam 2>{}/{}.isoseq3.polished.log".format(refParams.threads, subreadsXML, logDir, paramObj.group)
    subprocess.call(cmd, shell=True)

    cmd = "samtools fastq -@ {} polished.flnc.bam > polished.flnc.fq".format(refParams.threads)
    subprocess.call(cmd, shell=True)

    proovreadCorr(seqNeedCorrect="polished.flnc.fq", paramObj=paramObj, refParams=refParams, projectsDict=projectsDict)
    proovreadFq = os.path.join(os.getcwd(), "proovread", "proovread.fq")

    cmd = "minimap2 -ax splice:hq -uf {} -t {} --secondary=no --MD {} >aln.sam 2>{}/{}.polish.flnc.mm2.log".format(
        proovreadFq, refParams.ref_mm2_index, refParams.threads, logDir, paramObj.group)
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    if os.path.islink("aln.sam") or os.path.exists("aln.sam"):
        os.remove("aln.sam")
    mappedSam = os.path.join(prevDir, "isoseq3", "aln.sam")
    makeLink(mappedSam, "aln.sam")
    print str(datetime.datetime.now()) + " IsoSeq3 correction done!"

def fmlrc2Corr(tgsSample=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Fmlrc correction for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    resolveDir("fmlrc")
    logDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "log")
    leftReadsList = [i.strip() for i in re.split("[;|,]", tgsSample.ngsLeftReads)]
    rightReadsList = [i.strip() for i in re.split("[;|,]", tgsSample.ngsRightReads)]
    cmd = "zcat {} > leftReads.fastq".format(" ".join(leftReadsList))
    subprocess.call(cmd, shell=True)
    cmd = "zcat {} > rightReads.fastq".format(" ".join(rightReadsList))
    subprocess.call(cmd ,shell=True)
    cmd = "cat leftReads.fastq rightReads.fastq | awk 'NR % 4 == 2' | tr NT TN | "
    cmd = cmd + "ropebwt2 -LR 2>{}/{}.ropebwt2.log | tr NT TN | fmlrc2-convert comp_msbwt.npy 1>{}/{}.fmlrc_convert.log 2>&1".format(logDir, tgsSample.sampleName, logDir, tgsSample.sampleName)
    subprocess.call(cmd, shell=True)
    cmd = "seqkit seq --rna2dna {} -w 0 > nanopore.dna.fastq".format(tgsSample.tgsProcessedData)
    subprocess.call(cmd, shell=True)
    cmd = "fmlrc2 -t {} -C 10 comp_msbwt.npy nanopore.dna.fastq nanopore.fmlrc_corrected.fasta 1>{}/{}.fmlrc.log 2>&1".format(tgsSample.threads, logDir, tgsSample.sampleName)
    subprocess.call(cmd, shell=True)
    tgsSample.tgsProcessedData = os.path.join(os.getcwd(), "nanopore.fmlrc_corrected.fasta")
    # removeFiles(os.getcwd(), ["leftReads.fastq", "rightReads.fastq"])
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Fmlrc correction for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)


def mecatCorr(nanoporeRawSeq, tgsSample=None, refParams=None, mecatParams=None, useProovread=False, dirSpec=None):
    print str(datetime.datetime.now()) + " MECAT correction for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    # resolveDir("mecat")

    logDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "log")
    makeLink(nanoporeRawSeq, "nanoporeRawSeq.fq")
    nanoporeRawSeq = os.path.join(os.getcwd(), "nanoporeRawSeq.fq")
    if useProovread:
        proovreadCorr(nanoporeRawSeq, tgsSample=tgsSample, dirSpec=dirSpec)
        nanoporeRawSeq = os.path.join(os.getcwd(), "proovread", "proovread.fq")
    refGenomeSize = os.path.getsize(refParams.ref_genome)
    if refGenomeSize >= 200:
        mecatChunkSize = "50"
    elif refGenomeSize < 200 and refGenomeSize > 20:
        mecatChunkSize = "100"
    else:
        mecatChunkSize = "200"

    cmd = "mecat2pw -j 0 -d {} -o condidate.txt -w wrk_dir" \
          " -t {} -x {} -n {} -a {}".format(nanoporeRawSeq, tgsSample.threads, mecatParams.mecat_mode, mecatChunkSize,
                                            mecatParams.min_overlap)
    subprocess.call(cmd, shell=True)
    cmd = "mecat2cns -i 0 -t {} -x {} -p {} -r {} -a {} -c {} -l {} condidate.txt {} corrected.fa".format(
        tgsSample.threads, mecatParams.mecat_mode, mecatParams.part_num, mecatParams.min_read_score,
        mecatParams.min_overlap, mecatParams.min_read_cov, mecatParams.min_read_length, nanoporeRawSeq)
    subprocess.call(cmd, shell=True)
    tgsSample.tgsProcessedData = os.path.join(prevDir, "corrected.fa")

    # cmd = "minimap2 -ax splice -G 150k -k14 -uf {} --secondary=no corrected.fa -t {} >aln.mecat.sam 2>{}/{}.mecat.ont.mm2.log".format(
    #     refParams.ref_mm2_index, tgsSample.threads, logDir, tgsSample.sampleName)
    # subprocess.call(cmd, shell=True)
    # os.chdir(prevDir)
    # if os.path.islink("aln.sam") or os.path.exists("aln.sam"):
    #     os.remove("aln.sam")
    # mappedSam = os.path.join(prevDir, "mecat", "aln.sam")
    # makeLink(mappedSam, "aln.sam")
    print str(datetime.datetime.now()) + " MECAT correction for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)

def processReadsFromPacbio(paramObj=None, refParams=None, ccsParams=None):
    tgs_dir, strategy, projectName = paramObj.tgsDataDir, paramObj.tgsStrategy, paramObj.projectName
    # 1. Extract reads
    # 1.1 Determine the sequencing strategy
    logDir = os.path.join(refParams.out_dir, paramObj.projectName, paramObj.group, "log")
    if not os.path.exists(logDir):
        os.makedirs(logDir)
    # primer = determinePacPrimer(paramObj.tgsDataDir)
    print str(datetime.datetime.now()) + " Extract reads for project {} group {}...".format(paramObj.projectName, paramObj.group)
    if strategy.lower() == "rsii":
        print str(datetime.datetime.now()) + " The reads in {} are from RSII".format(paramObj.group)
        cmd = "find {} -name '*.bax.h5' > input.fofn".format(tgs_dir)
        subprocess.call(cmd, shell=True)
        cmd = "bax2bam -f input.fofn -o input"
        subprocess.call(cmd, shell=True)

        inputBam = "input.subreads.bam"
        print str(datetime.datetime.now()) + " Extract CCS reads from bax.h5 files..."
        cmd = "ccs {} CCS.bam -j {} --noPolish --minPasses {} --minLength {} --minPredictedAccuracy {} --minReadScore " \
              "{}".format(inputBam, refParams.threads, ccsParams.min_pass, ccsParams.min_subread_length,
                          ccsParams.min_predicted_accuracy, ccsParams.min_read_score)
        subprocess.call(cmd, shell=True)
        print str(datetime.datetime.now()) + " Extract CCS reads from bax.h5 files for project {} group {} done!".format(paramObj.projectName, paramObj.group)
    elif strategy.lower() == "sequel":
        print str(datetime.datetime.now()) + " The reads in project {} group {} are from Sequel".format(paramObj.projectName, paramObj.group)
        xmlFiles = os.popen("find -L {} -name '*.subreadset.xml'".format(tgs_dir))
        xmlFileStr = " ".join([i.strip("\n") for i in xmlFiles])
        if not xmlFileStr:
            bamFiles = os.popen("find -L {} -name '*.bam'".format(tgs_dir))
            bamFileStr = " ".join([i.strip("\n") for i in bamFiles])
            if not bamFileStr:
                raise Exception("You should put the *.bam in project {} the {} directory".format(paramObj.projectName, tgs_dir))
            fileStr = bamFileStr
        else:
            fileStr = xmlFileStr
        subreads = "input.subreadset.xml"
        cmd = "dataset create --type SubreadSet --force --name {} {} {}".format(projectName, subreads, fileStr)
        subprocess.call(cmd, shell=True)
        print str(datetime.datetime.now()) + " Extract CCS reads from subreads files for project {} group {}...".format(paramObj.projectName, paramObj.group)
        ccsOut = subprocess.Popen(["ccs", "--version"], stdout=subprocess.PIPE).communicate()[0].strip("\n")
        if "ccs 3." in ccsOut:
            cmd = "ccs {} CCS.bam -j {} --minPasses {} --minLength {} --minPredictedAccuracy {} --minReadScore " \
                  "{} --force".format(subreads, refParams.threads, ccsParams.min_pass, ccsParams.min_subread_length,
                              ccsParams.min_predicted_accuracy, ccsParams.min_read_score)
            subprocess.call(cmd, shell=True)
        elif "ccs 4." in ccsOut:
            cmd = "ccs {} CCS.bam -j {} --min-passes {} --min-length {} --min-rq {}".format(subreads, refParams.threads, ccsParams.min_pass, ccsParams.min_subread_length, ccsParams.min_predicted_accuracy)
            subprocess.call(cmd, shell=True)
        else:
            raise Exception("You should use version of ccs up to 3.4 or higher")
        print str(datetime.datetime.now()) + " Extract CCS reads from subreads files for project {} group {} done!".format(paramObj.projectName, paramObj.group)
    else:
        raise Exception("You should input the right PacBio sequencing strategy!!")

    inputBam = "CCS.bam"
    return inputBam

def processReadsFromNGS(tgsSample=None, refParams=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Start RNA-seq reads processing for project {} entry {}...".format(
        tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    # projectName = paramObj.projectName
    # groupName = paramObj.group
    resolveDir(os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "RNA-seq"))
    hisat2indexDir = os.path.join(dirSpec.out_dir, "hisat2indexDir")
    batchMakeDir(["alignment", "reassembly", "quant"])

    refAnnoGTF = refParams.ref_gtf
    gtfPrefix = os.path.splitext(os.path.basename(refParams.ref_gtf))[0]

    logDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "log")
    os.chdir("alignment")

    mergeList = []
    bamList = []
    juncList = []
    hisat2mappingLogList = []
    if tgsSample.ngsReadPair == "paired":
        leftReadsRepeats = tgsSample.ngsLeftReads.split(";")
        rightReadsRepeats = tgsSample.ngsRightReads.split(";")
        for i in range(len(leftReadsRepeats)):
            leftReads = ",".join([r.strip() for r in leftReadsRepeats[i].split(",")])
            rightReads = ",".join([r.strip() for r in rightReadsRepeats[i].split(",")])
            repeatName = "repeat" + str(i)
            resolveDir(repeatName)
            cmd = "hisat2 -x {}/{} -1 {} -2 {} --dta -p {} --max-intronlen 150000 --novel-splicesite-outfile {}.ss " \
              "--un-conc-gz {}.unmapped.fastq.gz 2>{}/{}.{}.hisat2.log | samtools sort -@ {} -o {}.sorted.bam ".format(
            hisat2indexDir, gtfPrefix, leftReads, rightReads, tgsSample.threads, repeatName, repeatName, logDir,
            tgsSample.sampleName, repeatName, tgsSample.threads, repeatName)
            subprocess.call(cmd, shell=True)

            cmd = "stringtie {}.sorted.bam -o {}.gtf -p {}".format(repeatName, repeatName, tgsSample.threads)
            subprocess.call(cmd, shell=True)
            # pysam.index("{}.sorted.bam".format(n.sampleName), catch_stdout=False)
            cmd = "samtools index -@ {} {}.sorted.bam".format(tgsSample.threads, repeatName)
            subprocess.call(cmd, shell=True)
            bamList.append("{}/{}.sorted.bam".format(os.getcwd(), repeatName))
            mergeList.append("{}/{}.gtf".format(os.getcwd(), repeatName))
            juncList.append("{}/{}.ss".format(os.getcwd(), repeatName))
            hisat2mappingLogList.append("{}/{}.{}.hisat2.log".format(logDir, tgsSample.sampleName, repeatName))
            os.chdir("../")
    else:
        if tgsSample.leftReads and tgsSample.rightReads == None:
            singleReadsRepeats = tgsSample.leftReads.split(";")
        elif tgsSample.leftReads == None and tgsSample.rightReads:
            singleReadsRepeats = tgsSample.rightReads.split(";")
        else:
            raise Exception("The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

        for i in range(len(singleReadsRepeats)):
            singleReads = ",".join([i.strip() for i in singleReadsRepeats[i].split(",")])
            repeatName = "repeat" + str(i)
            resolveDir(repeatName)
            cmd = "hisat2 -x {}/{} -U {} --dta -p {} --max-intronlen 150000 --novel-splicesite-outfile {}.ss " \
                  "--un-conc {}.unmmaped.fastq 2>{}/{}.{}.hisat2.log | samtools sort -@ {} -o {}.sorted.bam".format(
                hisat2indexDir, gtfPrefix, singleReads, tgsSample.threads, repeatName, repeatName, logDir,
                tgsSample.sampleName, repeatName, tgsSample.threads, repeatName
            )
            subprocess.call(cmd, shell=True)
            cmd = "stringtie {}.sorted.bam -o {}.gtf -p {}".format(repeatName, repeatName, tgsSample.threads)
            subprocess.call(cmd, shell=True)
            cmd = "samtools index -@ {} {}.sorted.bam".format(tgsSample.threads, repeatName)
            subprocess.call(cmd, shell=True)
            bamList.append("{}/{}.sorted.bam".format(os.getcwd(), repeatName))
            mergeList.append("{}/{}.gtf".format(os.getcwd(), repeatName))
            juncList.append("{}/{}.ss".format(os.getcwd(), repeatName))
            hisat2mappingLogList.append("{}/{}.{}.hisat2.log".format(logDir, tgsSample.sampleName, repeatName))
            os.chdir("../")

    os.chdir("../quant")
    bamStr = " ".join(bamList)
    if tgsSample.threads > 64:
        featureCountThreads = 64
    else:
        featureCountThreads = tgsSample.threads
    cmd = "featureCounts -a {} -T {} -o featureCounts_output.txt -G {} -J {} 1>{}/{}.featureCounts.log 2>&1".format(
        refAnnoGTF, featureCountThreads, refParams.ref_genome, bamStr, logDir, tgsSample.sampleName)
    subprocess.call(cmd, shell=True)
    os.chdir("../reassembly")
    gtfStr = " ".join(mergeList)
    cmd = "stringtie --merge -p {} -o stringtie_merged.gtf {}".format(tgsSample.threads, gtfStr)
    subprocess.call(cmd, shell=True)
    cmd = "samtools cat {} | samtools sort -@ {} > tmp.bam && samtools index tmp.bam".format(bamStr, tgsSample.threads)
    subprocess.call(cmd, shell=True)
    allMappedReadCounts = sum([getCountsFromHisat2MappingLog(i) for i in hisat2mappingLogList])
    getJuncCountFromNGSassembly("stringtie_merged.gtf", "tmp.bam", allReadCounts=allMappedReadCounts,
                                outFile="junctions.bed12", threads=tgsSample.threads, filterCount=3, juncList=juncList)
    tgsSample.ngsJunctions = os.path.join(os.getcwd(), "junctions.bed12")
    # removeFiles(os.getcwd(), ["tmp.bam", "tmp.bam.bai"])
    print str(datetime.datetime.now()) + " RNA-seq reads processing for project {} entry {} done!".format(
        tgsSample.projectName, tgsSample.sampleName)
    os.chdir(prevDir)

def proovreadCorr(seqNeedCorrect, tgsSample=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Proovread correction for project {} entry {}...".format(paramObj.projectName, paramObj.sampleName)
    prevDir = os.getcwd()
    resolveDir("proovread")
    logDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "log")
    if not os.path.exists(logDir):
        os.makedirs(logDir)
    if not os.path.exists("tmp"):
        os.makedirs("tmp")

    projectName, sampleName = tgsSample.projectName, tgsSample.sampleName

    maxChunkSize = str(tgsSample.memory[0:-1]) / 20
    if maxChunkSize <= 1000:
        proovreadChunkSize = maxChunkSize
    else:
        proovreadChunkSize = "1G"

    cmd = "SeqChunker -s {} -o pb-%03d.fq {}".format(proovreadChunkSize, seqNeedCorrect)
    subprocess.call(cmd, shell=True)
    fqFiles = os.popen("ls pb-*.fq")

    leftReadsList = [i.strip() for i in re.split("[;|,]", tgsSample.ngsLeftReads)]
    rightReadsList = [i.strip() for i in re.split("[;|,]", tgsSample.ngsRightReads)]
    cmd = "gunzip -k -c {} > leftReads.fastq".format(" ".join(leftReadsList))
    subprocess.call(cmd, shell=True)
    cmd = "gunzip -k -c {} > rightReads.fastq".format(" ".join(rightReadsList))
    subprocess.call(cmd, shell=True)
    leftReads, rightReads = "leftReads.fastq", "leftReads.fastq"
    for i in fqFiles:
        chunkName = i.strip("\n").split(".")[0]
        cmd = "proovread -l {} -s {} -s {} --prefix {} -t {} --overwrite " \
              "2>{}/{}.{}.proovread.log".format(i.strip("\n"), leftReads, rightReads, chunkName, tgsSample.threads, logDir, sampleName, chunkName)
        subprocess.call(cmd, shell=True)
        os.chdir(chunkName)
        cmd = '''seqkit grep ../{} -f <(cut -f 1 {}.ignored.tsv) -w 0 | cat {}.untrimmed.fq - > {}.proovread.fq'''.format(i.strip("\n"), chunkName, chunkName, chunkName)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        os.chdir("../")
        cmd = "cd tmp/ && ln -sf ../{}/{}.proovread.fq ./ && cd ../".format(chunkName, chunkName)
        subprocess.call(cmd, shell=True)
    cmd = "rm pb-*.fq && rm leftReads.fastq rightReads.fastq"
    # subprocess.call(cmd, shell=True)

    cmd = "cat tmp/*.proovread.fq > proovread.fq"
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Proovread correction for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)
