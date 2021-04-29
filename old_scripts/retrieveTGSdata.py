#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: retrieveTGSdata.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-11-27 22:16:58
Last modified: 2020-11-27 22:16:59
'''

from commonFuncs import *

def processReadsFromNanopore(tmpRun, refParams, threads, workDir):
    resolveDir(workDir)
    tgsDir = tmpRun.dataLocation
    nanoporeFileType = checkFast5Files(tgsDir)
    if nanoporeFileType == "fast5":
        flowcellType = "FLO-MIN106"
        kitType = "SQK-RNA002"
        basecallingOut = "albacore_basecalling"
        cmd = "read_fast5_basecaller.py --input {} --recursive --worker_threads {} --flowcell {} --kit {}" \
              " --save_path {}".format(tgsDir, threads, flowcellType, kitType, basecallingOut)
        subprocess.call(cmd, shell=True)
        nanoporeRawFq = "nanoporeRawSeq.fq"
        cmd = "cat {}/workspace/pass/* > {}".format(basecallingOut, nanoporeRawFq)
        subprocess.call(cmd, shell=True)
        cmd = "nanopolish index --directory={} --sequencing-summary={}/sequencing_summary.txt {}".format(tgsDir, basecallingOut, nanoporeRawFq)
        subprocess.call(cmd, shell=True)
        if not os.path.exists(refParams.ref_mm2_index):
            cmd = "minimap2 -d {} {}".format(refParams.ref_mm2_index, refParams.ref_genome)
            subprocess.call(cmd, shell=True)
        cmd = "minimap2 -ax splice -G 150K -uf -k 14 -t {} --secondary=no {} {} | samtools view -b - -o nanopore.mapped.bam".format(tmpRun.threads, refParams.ref_mm2_index, nanoporeRawFq)
        subprocess.call(cmd, shell=True)
        cmd = "samtools sort nanopore.mapped.bam > nanopore.mapped.sorted.bam && samtools index nanopore.mapped.sort.bam"
        subprocess.call(cmd, shell=True)
        cmd = "nanopolish polya --threads={} --reads={} --bam=nanopore.mapped.sorted.bam --genome={} > polya_results.tsv".format(threads, nanoporeRawFq, refParams.ref_genome)
        subprocess.call(cmd, shell=True)
        return os.path.join(os.getcwd(), nanoporeRawFq)
    else:
        cmd = "cat {}/* > nanoporeRawSeq.fq".format(tgsDir)
        subprocess.call(cmd, shell=True)
        return os.path.join(os.getcwd(), "nanoporeRawSeq.fq")

def retrieveTGSdata(tgsPlat2sample, dirSpec, refInfoParams, ccsParams, optionTools):
    projectNum = sum([len(tgsPlat2sample[i][j]) for i in tgsPlat2sample for j in tgsPlat2sample[i]])
    threads = optionTools.threads / projectNum

    pool = MyPool(processes=projectNum)
    for projectName in tgsPlat2sample:
        retrieveDir = os.path.join(dirSpec.out_dir, projectName, "retrieveData")
        for tgsPlat in tgsPlat2sample[projectName]:
            if tgsPlat == "pacbio":
                for tmpRun in tgsPlat2sample[projectName][tgsPlat]:
                    pacbioDemultiplexDir = os.path.join(retrieveDir, "pacbioDemultiplexed")
                    # runCCSandLima(tmpRun, ccsParams, threads, pacbioDemultiplexDir)
                    resolveDir(pacbioDemultiplexDir, chdir=False)
                    pool.apply_async(runCCSandLima, (tmpRun, ccsParams, threads, pacbioDemultiplexDir))
            elif tgsPlat == "nanopore":
                for tmpRun in tgsPlat2sample[projectName][tgsPlat]:
                    workDir = os.path.join(retrieveDir, "nanoporeRetrieve", tmpRun.strainName + "_" + tmpRun.groupedSampleName)
                    refParams = refInfoParams[tmpRun.refStrain]
                    pool.apply_async(processReadsFromNanopore, (tmpRun, refParams, threads, workDir))
            else:
                raise Exception("The TGS platform you input can't be identified, please check it!")
    pool.close()
    pool.join()

    for projectName in tgsPlat2sample:
        for tgsPlat in tgsPlat2sample[projectName]:
            if tgsPlat == "pacbio":
                demultiplexDir = os.path.join(dirSpec.out_dir, projectName, "retrieveData", "pacbioDemultiplexed")
                for sample in tgsPlat2sample[projectName][tgsPlat]:
                    groupedSampleName = sample.groupedSampleName
                    primers = sample.primers
                    records = SeqIO.parse(primers, "fasta")
                    recordsDict = SeqIO.to_dict(records)
                    for i in groupedSampleName.split("+"):
                        primerName = i
                        # tmpBamName = "{}_{}.fl.5p--{}_3p.bam".format(sample.projectName, sample.uniqName, primerName)
                        # cmd = "ls {}/{} > {}/{}.all.fofn".format(demultiplexDir, tmpBamName, demultiplexDir, sample.uniqName)
                        # subprocess.call(cmd, shell=True)

                        ids = ["5p", primerName + "_3p"]
                        resultRecords = [recordsDict[id_] for id_ in ids]
                        primerFileName = "{}.{}.primers.fa".format(sample.uniqName, primerName)
                        primerOut = open(os.path.join(demultiplexDir, primerFileName), "w")
                        SeqIO.write(resultRecords, primerOut, "fasta")
                        primerOut.close()

    pool = MyPool(processes=projectNum)
    for projectName in tgsPlat2sample:
        for tgsPlat in tgsPlat2sample[projectName]:
            if tgsPlat == "pacbio":
                demultiplexDir = os.path.join(dirSpec.out_dir, projectName, "retrieveData", "pacbioDemultiplexed")
                for sample in tgsPlat2sample[projectName][tgsPlat]:
                    groupedSampleName = sample.groupedSampleName
                    for i in groupedSampleName.split("+"):
                        primerName = i
                        # runRefine(primerName, demultiplexDir, threads)
                        # tmpBamName = "{}.fl.5p--{}_3p.bam".format(sample.uniqName, primerName)
                        # primerFileName = "{}.{}.primers.fa".format(sample.uniqName, primerName)
                        pool.apply_async(runRefine, (sample, primerName, demultiplexDir, threads))
    pool.close()
    pool.join()

def runCCSandLima(tmpRun, ccsParams, threads, workDir):
    projectName, strategy = tmpRun.projectName, tmpRun.tgsStrategy
    # tgs_dir, strategy, projectName = tmpRun.tgsDataDir, tmpRun.tgsStrategy, tmpRun.projectName
    # 1. Extract reads
    # 1.1 Determine the sequencing strategy
    print str(datetime.datetime.now()) + " Extract reads for project {} entry {}...".format(tmpRun.projectName,
                                                                                            tmpRun.groupedSampleName)
    resolveDir(workDir, chdir=True)

    if strategy.lower() == "rsii":
        print str(datetime.datetime.now()) + " The reads in project {} entry {} are from RSII".format(tmpRun.projectName, tmpRun.groupedSampleName)
        if not os.path.isdir(tmpRun.dataLocation):
            raise ValueError("The RSii data should be stored in a directory, please check it!")

        tgsDir = tmpRun.dataLocation
        cmd = "find {} -name '*.bax.h5' > {}.fofn".format(tgsDir, tmpRun.uniqName)
        subprocess.call(cmd, shell=True)
        cmd = "bax2bam -f {}.fofn -o {}".format(tmpRun.uniqName, tmpRun.uniqName)
        subprocess.call(cmd, shell=True)

        inputBam = "{}.subreads.bam".format(tmpRun.uniqName)
        print str(datetime.datetime.now()) + " Extract CCS reads from bax.h5 files and demultiplexing for project {} entry {}...".format(
            tmpRun.projectName, tmpRun.groupedSampleName)
        ccsOut = subprocess.Popen(["ccs", "--version"], stdout=subprocess.PIPE).communicate()[0].strip("\n")
        if "ccs 3." not in ccsOut:
            raise ValueError("The RSii data only support pbccs version 3.*, please use the correct version")

        cmd = "ccs {} {}.CCS.bam -j {} --noPolish --minPasses {} --minLength {} --minPredictedAccuracy {} --minReadScore " \
              "{}".format(inputBam, tmpRun.uniqName, threads, ccsParams.min_pass, ccsParams.min_subread_length,
                          ccsParams.min_predicted_accuracy, ccsParams.min_read_score)
        subprocess.call(cmd, shell=True)
        print str(datetime.datetime.now()) + " Extract CCS reads from bax.h5 files for project {} entry {} done!".format(
            tmpRun.projectName, tmpRun.groupedSampleName)
    elif strategy.lower() == "sequel":
        print str(datetime.datetime.now()) + " The reads in project {} entry {} are from Sequel".format(
            tmpRun.projectName, tmpRun.groupedSampleName)
        if os.path.isdir(tmpRun.dataLocation):
            tgsDir = tmpRun.dataLocation
            xmlFiles = os.popen("find -L {} -name '*.subreadset.xml'".format(tgsDir))
            xmlFileStr = " ".join([i.strip("\n") for i in xmlFiles])
            if not xmlFileStr:
                bamFiles = os.popen("find -L {} -name '*.bam'".format(tgsDir))
                bamFileStr = " ".join([i.strip("\n") for i in bamFiles])
                if not bamFileStr:
                    raise Exception(
                        "You should put the *.bam in project {} the {} directory".format(tmpRun.projectName, tgsDir))
                fileStr = bamFileStr
            else:
                fileStr = xmlFileStr
        else:
            tgsFiles = [i.strip() for i in tmpRun.dataLocation.split(",")]
            for i in tgsFiles:
                if not os.path.exists(i+".pbi"):
                    cmd = "pbindex {}".format(i)
                    subprocess.call(cmd, shell=True)
            tgsFilesStr = " ".join(tgsFiles)
            fileStr = tgsFilesStr

        subreads = "{}.subreadset.xml".format(tmpRun.uniqName)
        cmd = "dataset create --type SubreadSet --force --name {} {} {}".format(projectName, subreads, fileStr)
        subprocess.call(cmd, shell=True)

        print str(datetime.datetime.now()) + " Extract CCS reads from subreads files for project {} entry {}...".format(
            tmpRun.projectName, tmpRun.groupedSampleName)
        ccsOut = subprocess.Popen(["ccs", "--version"], stdout=subprocess.PIPE).communicate()[0].strip("\n")
        if "ccs 3." in ccsOut:
            cmd = "ccs {} {}.CCS.bam -j {} --minPasses {} --minLength {} --minPredictedAccuracy {} --minReadScore " \
                  "{} --force".format(subreads, tmpRun.uniqName, threads, ccsParams.min_pass, ccsParams.min_subread_length,
                                      ccsParams.min_predicted_accuracy, ccsParams.min_read_score)
            subprocess.call(cmd, shell=True)
        elif "ccs 4." in ccsOut:
            cmd = "ccs {} {}.CCS.bam -j {} --min-passes {} --min-length {} --min-rq {}".format(subreads, tmpRun.uniqName,
                                                                                               threads, ccsParams.min_pass,
                                                                                               ccsParams.min_subread_length,
                                                                                               ccsParams.min_predicted_accuracy)
            subprocess.call(cmd, shell=True)
        else:
            raise Exception("You should use version of ccs up to 3.4 or higher")
        print str(
            datetime.datetime.now()) + " Extract CCS reads from subreads files for project {} entry {} done!".format(
            tmpRun.projectName, tmpRun.groupedSampleName)
    else:
        raise Exception("You should input the right PacBio sequencing strategy!!")

    print str(datetime.datetime.now()) + " Demultiplex CCS bam for project {} entry {}...".format(tmpRun.projectName, tmpRun.groupedSampleName)
    cmd = "lima {}.CCS.bam {} {}.fl.bam --isoseq --dump-clips -j {}".format(tmpRun.uniqName, tmpRun.primersUsed, tmpRun.uniqName, threads)
    subprocess.call(cmd, shell=True)
    print str(datetime.datetime.now()) + " Demultiplex CCS bam for project {} entry {} done!".format(tmpRun.projectName,tmpRun.groupedSampleName)
    # inputBam = "CCS.bam"
    # return inputBam

def runRefine(tgsSample, primerName, demultiplexDir, threads):
    resolveDir(demultiplexDir)
    tmpBamName = "{}.fl.5p--{}_3p.bam".format(tgsSample.uniqName, primerName)
    primerFileName = "{}.{}.primers.fa".format(tgsSample.uniqName, primerName)
    cmd = "isoseq3 refine {} {} {}.{}.flnc.bam --require-polya -j {}".format(tmpBamName, primerFileName, tgsSample.uniqName, primerName, threads)
    subprocess.call(cmd, shell=True)
