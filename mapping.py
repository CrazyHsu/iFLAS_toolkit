#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: preprocess.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:11:47
Last modified: 2021-04-29 16:11:47
'''

from commonFuncs import *
from commonObjs import *

def getJuncFromRegtools(dataObj=None, dirSpec=None, filterByCount=10):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Get junctions from RNA-seq with Regtools for project {} sample {}...".format(projectName, sampleName)
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    logDir = os.path.join(baseDir, "log")
    resolveDir(logDir, chdir=False)
    rnaseqSortedBam = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "tmp.bam")
    cmd = "regtools junctions extract -a 5 -m 50 -M 50000 {} -s 0 -o tmp.bed 2>{}/regtools.log".format(rnaseqSortedBam, logDir)
    subprocess.call(cmd, shell=True)
    cmd = '''awk '{if($5>%d){print}}' tmp.bed > junctions.bed''' % (filterByCount)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    dataObj.ngs_junctions = os.path.join(os.getcwd(), "junctions.bed")
    removeFiles(os.getcwd(), ["tmp.bed"])
    print getCurrentTime() + " Get junctions from RNA-seq with Regtools for project {} sample {} done!".format(projectName, sampleName)

def filterByJunc(bedFile=None, refBedFile=None, juncBedFile=None, maxIntronLen=50000):
    refBed = BedFile(refBedFile, type="bed12+").reads
    flncBed = BedFile(bedFile, type="bed12+").reads

    ngsJuncDict = {}
    if juncBedFile:
        ngsJuncBed = BedFile(juncBedFile, type="bed12").reads
        for read in ngsJuncBed:
            if len(ngsJuncBed[read].introns) > 0:
                for junc in ngsJuncBed[read].introns:
                    if junc not in ngsJuncDict:
                        ngsJuncDict[junc] = ""
    refJuncDict = {}
    for trans in refBed:
        if len(refBed[trans].introns) > 0:
            for junc in refBed[trans].introns:
                if junc not in refJuncDict:
                    refJuncDict[junc] = ""

    flncJuncDict = {}
    for read in flncBed:
        if len(flncBed[read].introns) > 0:
            for junc in flncBed[read].introns:
                if junc not in flncJuncDict:
                    flncJuncDict[junc] = 1
                else:
                    flncJuncDict[junc] += 1

    filterReads = open("filterReads.lst", "w")
    for read in flncBed:
        if len(flncBed[read].introns) > 0:
            if getBlockLength(flncBed[read].introns) > maxIntronLen: continue
            consensusJuncN = 0
            for junc in flncBed[read].introns:
                if junc in refJuncDict or junc in ngsJuncDict or flncJuncDict[junc] >= 3:
                    consensusJuncN += 1
            if len(flncBed[read].introns) >= 4:
                print >> filterReads, read
            elif consensusJuncN/float(len(flncBed[read].introns)) >= 1.0/3:
                print >> filterReads, read
        else:
            print >> filterReads, read
    filterReads.close()
    return os.path.join(os.getcwd(), "filterReads.lst")

def mappingFilterAndAddTags(samFile=None, outPrefix="flnc", maxLength=50000, refBedFile=None, juncBedFile=None, threads=10):
    cmd = "bamToBed -i <(samtools view -bS {}) -bed12 > tmp.bed12".format(samFile)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    filteredReads = filterByJunc(bedFile="tmp.bed12", refBedFile=refBedFile, juncBedFile=juncBedFile, maxIntronLen=maxLength)
    cmd = '''(samtools view -H {}; filter.pl -o {} {} -m i) > tmp.sam'''.format(samFile, filteredReads, samFile)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''(samtools view -H tmp.sam; samtools view -f 16 -F 4079 tmp.sam; samtools view -f 0 -F 4095 tmp.sam) | 
        samAddTag.pl --checkHardClip --coverage --identity 2>lengthInconsistent.sam |
        samtools sort -@ {} --output-fmt SAM > {}.mm2.sam'''.format(threads, outPrefix)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "sam2bed.pl -t CV,ID {}.mm2.sam > {}.addCVandID.bed12+".format(outPrefix, outPrefix)
    subprocess.call(cmd, shell=True)
    cmd = "samtools view -@ {} -bS {}.mm2.sam > {}.mm2.sorted.bam".format(threads, outPrefix, outPrefix)
    subprocess.call(cmd, shell=True)
    removeFiles(os.getcwd(), ["tmp.bed12", "tmp.sam"])

def minimap2mapping(dataObj=None, minimap2Params=None, refParams=None, dirSpec=None, threads=10):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Mapping flnc reads to reference genome for project {} sample {}...".format(projectName, sampleName)
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    prevDir = os.getcwd()
    workDir = os.path.join(baseDir, "mapping")
    resolveDir(workDir)
    if minimap2Params.mm2_index != None and validateFile(minimap2Params.mm2_index):
        mm2index = minimap2Params.mm2_index
    elif refParams.ref_mm2_index != None and validateFile(refParams.ref_mm2_index):
        mm2index = refParams.ref_mm2_index
    else:
        cmd = "minimap2 -d {} -t {} {}".format("ref.mm2", threads, refParams.ref_genome)
        subprocess.call(cmd, shell=True)
        mm2index = os.path.join(workDir, "ref.mm2")
        dataObj.mm2index = mm2index
    sampleName = dataObj.sample_name
    processedFlncFq = dataObj.data_processed_location
    logDir = os.path.join(baseDir, "log")
    resolveDir(logDir, chdir=False)
    if dataObj.tgs_plat == "pacbio":
        cmd = "minimap2 -ax splice:hq -G {} -uf {} --secondary=no --MD {} -t {} >flnc.mm2.sam 2>{}/{}.mm2.log".format(
            minimap2Params.max_intron_length, mm2index, processedFlncFq, threads, logDir, sampleName)
        subprocess.call(cmd, shell=True)
        if dataObj.data_processed_location == None:
            rawFlncFq = os.path.join(baseDir, "preprocess", "pacbio", "rawFlnc.fq")
        else:
            rawFlncFq = os.path.join(baseDir, "preprocess", "fmlrc", "rawFlnc.fq")
        cmd = "minimap2 -ax splice:hq -G {} -uf {} --secondary=no --MD {} -t {} >rawFlnc.mm2.sam 2>{}/{}.rawFlnc.mm2.log".format(
            minimap2Params.max_intron_length, mm2index, rawFlncFq, threads, logDir, sampleName)
        subprocess.call(cmd, shell=True)
    else:
        cmd = "minimap2 -ax splice -G {} -k14 -uf {} --secondary=no --MD {} -t {} >flnc.mm2.sam 2>{}/{}.mm2.log".format(
            minimap2Params.max_intron_length, mm2index, processedFlncFq, threads, logDir, sampleName)
        subprocess.call(cmd, shell=True)
        if dataObj.data_processed_location == None:
            rawFlncFq = os.path.join(baseDir, "preprocess", "nanopore", "rawFlnc.fq")
        else:
            rawFlncFq = os.path.join(baseDir, "preprocess", "fmlrc", "rawFlnc.fq")
        cmd = "minimap2 -ax splice -G {} -k14 -uf {} --secondary=no --MD {} -t {} >rawFlnc.mm2.sam 2>{}/{}.rawFlnc.mm2.log".format(
            minimap2Params.max_intron_length, mm2index, rawFlncFq, threads, logDir, sampleName)
        subprocess.call(cmd, shell=True)

    cmd = "seqkit fq2fa -w 0 {} > flnc.processed.fa".format(processedFlncFq)
    subprocess.call(cmd, shell=True)
    if dataObj.ngs_junctions == None and (dataObj.ngs_right_reads or dataObj.ngs_left_reads):
        dataObj.ngs_junctions = os.path.join(baseDir, "mapping", "rna-seq", "reassembly", "junctions.bed")
    mappingFilterAndAddTags(samFile="flnc.mm2.sam", outPrefix="flnc", maxLength=minimap2Params.max_intron_length, refBedFile=refParams.ref_bed, juncBedFile=dataObj.ngs_junctions, threads=threads)
    mappingFilterAndAddTags(samFile="rawFlnc.mm2.sam", outPrefix="rawFlnc", maxLength=minimap2Params.max_intron_length, refBedFile=refParams.ref_bed, juncBedFile=dataObj.ngs_junctions, threads=threads)
    print getCurrentTime() + " Mapping flnc reads to reference genome for project {} sample {} done!".format(projectName, sampleName)
    os.chdir(prevDir)

def hisat2mapping(dataObj=None, refParams=None, dirSpec=None, threads=10):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Mapping rna-seq short reads to reference genome with hisat2 for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    workDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "mapping", "rna-seq")
    resolveDir(workDir)
    logDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "log")
    resolveDir(logDir, chdir=False)
    # hisat2indexDir = os.path.join(dirSpec.out_dir, "hisat2indexDir")
    # gtfPrefix = os.path.splitext(os.path.basename(refParams.ref_gtf))[0]
    checkAndMakeHisat2Index(refParams=refParams, dirSpec=dirSpec, threads=threads)

    batchCreateDir(["alignment", "reassembly"])
    resolveDir("alignment")
    bamList = []
    mergeList = []
    if dataObj.ngs_reads_paired == "paired":
        leftReadsRepeats = dataObj.ngs_left_reads.split(";")
        rightReadsRepeats = dataObj.ngs_right_reads.split(";")
        for i in range(len(leftReadsRepeats)):
            leftReads = ",".join([r.strip() for r in leftReadsRepeats[i].split(",")])
            rightReads = ",".join([r.strip() for r in rightReadsRepeats[i].split(",")])
            repeatName = "repeat" + str(i)
            resolveDir(repeatName)
            cmd = "hisat2 -x {} -1 {} -2 {} --dta -p {} --max-intronlen 50000 --novel-splicesite-outfile {}.ss " \
              "--un-conc-gz {}.unmapped.fastq.gz 2>{}/{}.{}.hisat2.log | samtools sort -@ {} -o {}.sorted.bam ".format(
            refParams.hisat2_index, leftReads, rightReads, threads, repeatName, repeatName, logDir,
            sampleName, repeatName, threads, repeatName)
            subprocess.call(cmd, shell=True)

            cmd = "stringtie {}.sorted.bam -o {}.gtf -p {}".format(repeatName, repeatName, threads)
            subprocess.call(cmd, shell=True)
            cmd = "samtools index -@ {} {}.sorted.bam".format(threads, repeatName)
            subprocess.call(cmd, shell=True)
            bamList.append("{}/{}.sorted.bam".format(os.getcwd(), repeatName))
            mergeList.append("{}/{}.gtf".format(os.getcwd(), repeatName))
            os.chdir("../")
    else:
        if dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
            singleReadsRepeats = dataObj.ngs_left_reads.split(";")
        elif dataObj.ngs_left_reads == None and dataObj.ngs_right_reads:
            singleReadsRepeats = dataObj.ngs_right_reads.split(";")
        else:
            raise Exception("The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

        for i in range(len(singleReadsRepeats)):
            singleReads = ",".join([i.strip() for i in singleReadsRepeats[i].split(",")])
            repeatName = "repeat" + str(i)
            resolveDir(repeatName)
            cmd = "hisat2 -x {} -U {} --dta -p {} --max-intronlen 50000 --novel-splicesite-outfile {}.ss " \
                  "--un-conc {}.unmmaped.fastq 2>{}/{}.{}.hisat2.log | samtools sort -@ {} -o {}.sorted.bam".format(
                refParams.hisat_index, singleReads, threads, repeatName, repeatName, logDir,
                sampleName, repeatName, threads, repeatName)
            subprocess.call(cmd, shell=True)
            cmd = "stringtie {}.sorted.bam -o {}.gtf -p {}".format(repeatName, repeatName, threads)
            subprocess.call(cmd, shell=True)
            cmd = "samtools index -@ {} {}.sorted.bam".format(threads, repeatName)
            subprocess.call(cmd, shell=True)
            bamList.append("{}/{}.sorted.bam".format(os.getcwd(), repeatName))
            mergeList.append("{}/{}.gtf".format(os.getcwd(), repeatName))
            os.chdir("../")

    os.chdir(os.path.join(workDir, "reassembly"))
    gtfStr = " ".join(mergeList)
    bamStr = " ".join(bamList)
    cmd = "stringtie --merge -p {} -o stringtie_merged.gtf {}".format(threads, gtfStr)
    subprocess.call(cmd, shell=True)
    cmd = "samtools cat {} | samtools sort -@ {} > tmp.bam && samtools index tmp.bam".format(bamStr, threads)
    subprocess.call(cmd, shell=True)
    if dataObj.ngs_junctions == None:
        getJuncFromRegtools(dataObj=dataObj, dirSpec=dirSpec)
    resolveDir(prevDir)
    print getCurrentTime() + " Mapping rna-seq short reads to reference genome with hisat2 for project {} sample {} done!".format(projectName, sampleName)

def mapping(dataObj=None, minimap2Params=None, refParams=None, dirSpec=None, threads=10, useFmlrc2=True):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
        from preprocess import renameNGSdata2fastp, processRnaseq
        processRnaseq(dataObj=dataObj, threads=threads, dirSpec=dirSpec, max_reads_length_tirmmed=1)
        renameNGSdata2fastp(dataObj=dataObj)
        hisat2mapping(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, threads=threads)

    if dataObj.data_processed_location:
        if isinstance(dataObj.data_processed_location, list):
            validFiles = []
            for i in dataObj.data_processed_location:
                if validateFile(i):
                    validFiles.append(i)
            preprocessDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "preprocess", dataObj.tgs_plat.lower())
            resolveDir(preprocessDir, chdir=False)
            dataObj.data_processed_location = os.path.join(preprocessDir, "rawFlnc.fq")
            cmd = "cat {} > {}".format(" ".join(validFiles), dataObj.data_processed_location)
            subprocess.call(cmd, shell=True)
            if dataObj.use_fmlrc2 and useFmlrc2:
                from preprocess import correctWithFmlrc2
                correctWithFmlrc2(dataObj, dirSpec=dirSpec, useFmlrc2=True, threads=dataObj.single_run_threads)
            minimap2mapping(dataObj=dataObj, minimap2Params=minimap2Params, refParams=refParams, dirSpec=dirSpec,
                            threads=threads)
        elif isinstance(dataObj.data_processed_location, basestring):
            if validateFile(dataObj.data_processed_location):
                if dataObj.use_fmlrc2 and useFmlrc2:
                    preprocessDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "preprocess", dataObj.tgs_plat.lower())
                    resolveDir(preprocessDir, chdir=False)
                    makeLink(dataObj.data_processed_location, os.path.join(preprocessDir, "rawFlnc.fq"))
                    from preprocess import correctWithFmlrc2
                    correctWithFmlrc2(dataObj, dirSpec=dirSpec, useFmlrc2=True, threads=dataObj.single_run_threads)
                minimap2mapping(dataObj=dataObj, minimap2Params=minimap2Params, refParams=refParams, dirSpec=dirSpec,
                                threads=threads)
    else:
        if dataObj.use_fmlrc2:
            dataObj.data_processed_location = os.path.join(dirSpec.out_dir, projectName, sampleName, "preprocess", "fmlrc", "fmlrc_corrected.fasta")
        else:
            dataObj.data_processed_location = os.path.join(dirSpec.out_dir, projectName, sampleName, "preprocess", dataObj.tgs_plat.lower(), "rawFlnc.fq")
            if useFmlrc2:
                from preprocess import correctWithFmlrc2
                correctWithFmlrc2(dataObj, dirSpec=dirSpec, useFmlrc2=True, threads=dataObj.single_run_threads)
        if validateFile(dataObj.data_processed_location):
            minimap2mapping(dataObj=dataObj, minimap2Params=minimap2Params, refParams=refParams, dirSpec=dirSpec, threads=threads)
        else:
            raise Exception("Something wrong happened for generating preprocess flnc reads! Please check it!")

