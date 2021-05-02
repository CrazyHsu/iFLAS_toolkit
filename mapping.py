#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: preprocess.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:11:47
Last modified: 2021-04-29 16:11:47
'''

from commonFuncs import *

def minimap2mapping(dataObj=None, minimap2Params=None, refParams=None, dirSpec=None, threads=10):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Mapping flnc reads to reference genome for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    workDir = os.path.join(dirSpec.out_dir, "mapping")
    resolveDir(workDir)
    if dataObj.mm2index != None:
        mm2index = minimap2Params.mm2index
    elif refParams.ref_mm2_index != None:
        mm2index = refParams.ref_mm2_index
    else:
        cmd = "minimap2 -d {} -t {} {}".format("ref.mm2", threads, refParams.ref_genome)
        subprocess.call(cmd, shell=True)
        mm2index = os.path.join(workDir, "ref.mm2")
        dataObj.mm2index = mm2index
    sampleName = dataObj.sample_name
    processedFlncFq = dataObj.data_processed_location
    logDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "log")
    if dataObj.tgs_plat == "pacbio":
        cmd = "minimap2 -ax splice:hq -G {} -uf {} --secondary=no --MD {} -t {} >{}.mm2.sam 2>{}/{}.mm2.log".format(
            minimap2Params.max_reads_length, mm2index, processedFlncFq, threads, sampleName, logDir, sampleName)
    else:
        cmd = "minimap2 -ax splice -G {} -k14 -uf {} --secondary=no --MD {} -t {} >{}.mm2.sam 2>{}/{}.mm2.log".format(
            minimap2Params.max_reads_length, mm2index, processedFlncFq, threads, sampleName, logDir, sampleName)
    subprocess.call(cmd, shell=True)
    print getCurrentTime() + " Mapping flnc reads to reference genome for project {} sample {} done!".format(projectName, sampleName)
    os.chdir(prevDir)

def hisat2mapping(dataObj=None, refParams=None, dirSpec=None, optionTools=None, threads=10):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Mapping rna-seq short reads to reference genome with hisat2 for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    workDir = os.path.join(dirSpec.out_dir, "mapping", "rna-seq")
    resolveDir(workDir)
    logDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "log")
    hisat2indexDir = os.path.join(workDir, "hisat2indexDir")
    gtfPrefix = os.path.splitext(os.path.basename(refParams.ref_gtf))[0]
    makeHisat2Index(refParams=refParams, hisat2indexDir=hisat2indexDir, indexPrefix=gtfPrefix, optionTools=optionTools)

    batchCreateDir(["alignment", "reassembly"])
    resolveDir("alignment")
    bamList = []
    mergeList = []
    if dataObj.ngsReadPair == "paired":
        leftReadsRepeats = dataObj.ngsLeftReads.split(";")
        rightReadsRepeats = dataObj.ngsRightReads.split(";")
        for i in range(len(leftReadsRepeats)):
            leftReads = ",".join([r.strip() for r in leftReadsRepeats[i].split(",")])
            rightReads = ",".join([r.strip() for r in rightReadsRepeats[i].split(",")])
            repeatName = "repeat" + str(i)
            resolveDir(repeatName)
            cmd = "hisat2 -x {}/{} -1 {} -2 {} --dta -p {} --max-intronlen 50000 --novel-splicesite-outfile {}.ss " \
              "--un-conc-gz {}.unmapped.fastq.gz 2>{}/{}.{}.hisat2.log | samtools sort -@ {} -o {}.sorted.bam ".format(
            hisat2indexDir, gtfPrefix, leftReads, rightReads, threads, repeatName, repeatName, logDir,
            dataObj.sampleName, repeatName, threads, repeatName)
            subprocess.call(cmd, shell=True)

            cmd = "stringtie {}.sorted.bam -o {}.gtf -p {}".format(repeatName, repeatName, threads)
            subprocess.call(cmd, shell=True)
            cmd = "samtools index -@ {} {}.sorted.bam".format(threads, repeatName)
            subprocess.call(cmd, shell=True)
            bamList.append("{}/{}.sorted.bam".format(os.getcwd(), repeatName))
            mergeList.append("{}/{}.gtf".format(os.getcwd(), repeatName))
            os.chdir("../")
    else:
        if dataObj.leftReads and dataObj.rightReads == None:
            singleReadsRepeats = dataObj.leftReads.split(";")
        elif dataObj.leftReads == None and dataObj.rightReads:
            singleReadsRepeats = dataObj.rightReads.split(";")
        else:
            raise Exception("The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

        for i in range(len(singleReadsRepeats)):
            singleReads = ",".join([i.strip() for i in singleReadsRepeats[i].split(",")])
            repeatName = "repeat" + str(i)
            resolveDir(repeatName)
            cmd = "hisat2 -x {}/{} -U {} --dta -p {} --max-intronlen 50000 --novel-splicesite-outfile {}.ss " \
                  "--un-conc {}.unmmaped.fastq 2>{}/{}.{}.hisat2.log | samtools sort -@ {} -o {}.sorted.bam".format(
                hisat2indexDir, gtfPrefix, singleReads, threads, repeatName, repeatName, logDir,
                dataObj.sampleName, repeatName, threads, repeatName)
            subprocess.call(cmd, shell=True)
            cmd = "stringtie {}.sorted.bam -o {}.gtf -p {}".format(repeatName, repeatName, threads)
            subprocess.call(cmd, shell=True)
            cmd = "samtools index -@ {} {}.sorted.bam".format(threads, repeatName)
            subprocess.call(cmd, shell=True)
            bamList.append("{}/{}.sorted.bam".format(os.getcwd(), repeatName))
            mergeList.append("{}/{}.gtf".format(os.getcwd(), repeatName))
            os.chdir("../")

    os.chdir("../reassembly")
    gtfStr = " ".join(mergeList)
    bamStr = " ".join(bamList)
    cmd = "stringtie --merge -p {} -o stringtie_merged.gtf {}".format(threads, gtfStr)
    subprocess.call(cmd, shell=True)
    cmd = "samtools cat {} | samtools sort -@ {} > tmp.bam && samtools index tmp.bam".format(bamStr, threads)
    subprocess.call(cmd, shell=True)
    resolveDir(prevDir)
    print getCurrentTime() + " Mapping rna-seq short reads to reference genome with hisat2 for project {} sample {} done!".format(
        projectName, sampleName)

def mapping(dataObj=None, minimap2Params=None, refParams=None, dirSpec=None, threads=10):
    if dataObj.data_processed_location:
        minimap2mapping(dataObj=dataObj, minimap2Params=minimap2Params, refParams=refParams, dirSpec=dirSpec, threads=threads)
    if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
        hisat2mapping(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, threads=threads)
