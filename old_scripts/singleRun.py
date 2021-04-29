#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
File name: singleRun.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-08-29 15:37:15
Last modified: 2019-08-29 15:37:16
'''
import argparse, os, psutil, datetime, shutil, sys, random
from multiprocessing import cpu_count
from multiprocessing import Pool
import subprocess
import pysam
from Bio import SeqIO
from Config import *
from findAS import *
from commonObjs import *
######## Parser arguments #############
parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-p", dest="project", type=str, default="project",
                    help="The project name")
parser.add_argument("-g", dest="group", type=str, default="group",
                    help="The group which the analysis belongs to. For single running, this parameter can be missed")
parser.add_argument("-t", dest="tgs_dir", type=str, default="tgsDir",
                    help="The directory which cantain the TGS raw reads to be analysis")
parser.add_argument("-n", dest="ngs_dir", type=str, default=None,
                    help="The directory which contain the NGS raw reads to be analysis")
parser.add_argument("-plat", dest="tgs_plat", type=str, default="PacBio",
                    help="The sequence strategy to carry out TGS sequencing. Default: PacBio")
parser.add_argument("-pac_plat", dest="pac_plat", type=str, default="sequel",
                    help="Specify the strategy that used to carry out PacBio sequencing, Sequel or RSII. Default: Sequel")
parser.add_argument("-pair", dest="ngs_pair", action="store_true",
                    help="If NGS raw reads provided, specify this parameter to determine the reads of NGS is paired")
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
parser.add_argument("-cc", dest="command", action="store_true",
                    help="Use this parameter to specify whether the program runs from command line or from "
                         "config file. Note: when the sequence information specified from command line the program will "
                         "ignore the information in config file and only be executed once")
args = parser.parse_args()

thread = psutil.cpu_count()*3/4
mem = psutil.virtual_memory().total
# paramObj = ResearchSection("singleRun")
# defaultCfg = None
# for k, v in args._get_kwargs():
#     if k == "default_cfg":
#         defaultCfg = Config(args.default_cfg)
#     else:
#         paramObj.__dict__[k] = v
# print paramObj, defaultCfg
#
# refParams = defaultCfg.refParams
# mecatParams = defaultCfg.mecatParams
# ccsParams = defaultCfg.ccsParams
#
binPath = "/home/xufeng/xufeng/iso-seq/code/isoSeq.6.24"

def initSysSetting(refParams=None, paramObj=None):
    # subprocess.call("source activate isoseq3_pipeline")
    # subprocess.call(["ccs", "-h"])
    global originDir
    originDir = os.getcwd()
    if refParams.ref_genome == None:
        raise Exception("You must specify a genome file!")
    if refParams.ref_gtf == None:
        raise Exception("You must specify a gtf file!")
    refParams.ref_genome = os.path.join(originDir, refParams.ref_genome)
    refParams.ref_gtf = os.path.join(originDir, refParams.ref_gtf)
    validateFile(refParams.ref_genome)
    validateFile(refParams.ref_gtf)
    if refParams.ref_size == None:
        refParams.ref_size = os.path.join(originDir, "default.len.txt")
        sizeOut = open(refParams.ref_size)
        for seqRecord in SeqIO.parse(refParams.ref_genome, "fasta"):
            print >>sizeOut, "\t".join([seqRecord.id, seqRecord.seq])
        sizeOut.close()
    else:
        refParams.ref_size = os.path.join(originDir, refParams.ref_size)
    if refParams.ref_gpe == None:
        refParams.ref_gpe = os.path.join(originDir, "default.gpe")
        cmd = "gtfToGenePred {} {} -genePredExt".format(refParams.ref_gtf, refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
    else:
        refParams.ref_gpe = os.path.join(originDir, refParams.ref_gpe)
    if refParams.ref2bit == None:
        refParams.ref2bit = os.path.join(originDir, "default.ref2bit")
        cmd = "faToTwoBit {} {}".format(refParams.ref_genome, refParams.ref2bit)
        subprocess.call(cmd, shell=True)
    else:
        refParams.ref2bit = os.path.join(originDir, refParams.ref2bit)
    refParams.out_dir = os.path.join(originDir, refParams.out_dir)
    refParams.ref_mm2_index = os.path.join(originDir, refParams.ref_mm2_index)

    # global thread, mem, paramObj, refParams, mecatParams, ccsParams
    # thread = psutil.cpu_count()*3/4
    # mem = psutil.virtual_memory().total
    # paramObj = ResearchSection("singleRun")
    # defaultCfg = None
    # for k, v in args._get_kwargs():
    #     if k == "default_cfg":
    #         defaultCfg = Config(args.default_cfg)
    #     else:
    #         paramObj.__dict__[k] = v
    # print paramObj, defaultCfg
    #
    # refParams = defaultCfg.refParams
    # mecatParams = defaultCfg.mecatParams
    # ccsParams = defaultCfg.ccsParams

def proovreadCorr(seqNeedCorrect, paramObj=None, refParams=None, tgs_plat="pacbio"):
    print str(datetime.datetime.now()) + " Proovread correction..."
    prevDir = os.getcwd()
    os.makedirs("proovread/log")
    os.makedirs("proovread/tmp")
    os.chdir("proovread")

    maxChunkSize = mem / 1024 / 1024 / 20
    if maxChunkSize <= 1000:
        proovreadChunkSize = maxChunkSize * 4 / 5
        proovreadChunkSize = "{}M".format(proovreadChunkSize)
    else:
        proovreadChunkSize = "1G"

    cmd = "SeqChunker -s {} -l 1 -o pb-%03d.fq ../{}".format(proovreadChunkSize, seqNeedCorrect)
    subprocess.call(cmd, shell=True)
    fqFiles = os.popen("ls pb-*.fq")
    leftReads = "zcat {}/*.R1.fastq.gz | gunzip > leftReads.fastq.gz".format(paramObj.ngs_dir)
    rightReads = "zcat {}/*.R2.fastq.gz | gunzip > rightReads.fastq.gz".format(paramObj.ngs_dir)
    for i in fqFiles:
        chunkName = i.split(".")[0]
        cmd = "proovread -l {} -s {} -s {} --prefix {} -t {} " \
              "2>log/{}.proovread.log".format(i, leftReads, rightReads, chunkName, thread, chunkName)
        subprocess.call(cmd, shell=True)
        cmd = "cd tmp/ & ln -sf ../{}/*.trimmed.fa ./ & cd ../".format(chunkName)
        subprocess.call(cmd, shell=True)
    cmd = "rm pb-*.fq && rm leftReads.fastq.gz rightReads.fastq.gz"
    subprocess.call(cmd, shell=True)

    cmd = "cat tmp/*.trimmed.fa > proovread.fa"
    subprocess.call(cmd, shell=True)
    if tgs_plat == "pacbio":
        cmd = "minimap2 -ax splice:hq -uf {} -t {} proovread.fa > aln.sam".format(refParams.ref_mm2_index, thread)
    else:
        cmd = "minimap2 -ax splice -uf -k14 {} -t {} proovread.fa > aln.sam".format(refParams.ref_mm2_index, thread)
    subprocess.call(cmd, shell=True)

    os.chdir(prevDir)
    os.symlink(os.path.join(prevDir, "proovread", "aln.sam"), "aln.sam")
    print str(datetime.datetime.now()) + " Proovread correction done!"

def makeRefTranscriptome():
    pass

def mecatCorr(seqNeedCorrect, refParams=None, mecatParams=None, tgs_plat="nanopore", hybrid=False):
    print str(datetime.datetime.now()) + " MECAT correction..."
    prevDir = os.getcwd()
    os.makedirs("mecat")
    os.chdir("mecat")

    refGenomeSize = os.path.getsize(refParams.ref_genome)
    if refGenomeSize >= 200:
        mecatChunkSize = "50"
    elif refGenomeSize < 200 and refGenomeSize > 20:
        mecatChunkSize = "100"
    else:
        mecatChunkSize = "200"

    mecatMode = 1 if tgs_plat == "nanopore" else 0
    cmd = "mecat2pw -j 0 -d {} -o condidate.txt -w wrk_dir" \
          " -t {} -x {} -n {} -a {}".format(seqNeedCorrect, thread, mecatMode, mecatChunkSize, mecatParams.min_overlap)
    subprocess.call(cmd, shell=True)
    cmd = "mecat2cns -i 0 -t {} -x {} -p {} -r {} -a {} -c {} -l {} condidate.txt {} corrected.fa".format(
        thread, mecatParams.mecat_mode, mecatParams.part_num, mecatParams.min_read_score, mecatParams.min_overlap,
        mecatParams.min_read_cov, mecatParams.min_read_length, seqNeedCorrect)
    subprocess.call(cmd, shell=True)
    if hybrid:
        print str(datetime.datetime.now()) + " MECAT correction done!"
        return os.path.join(prevDir, "mecat", "corrected.fa")
    else:
        cmd = "minimap2 -ax splice -k14 -uf {} corrected.fa -t {} > aln.sam".format(refParams.ref_mm2_index, thread)
        subprocess.call(cmd, shell=True)
        os.chdir(prevDir)
        os.symlink(os.path.join(prevDir, "mecat", "aln.sam"), "aln.sam")
        print str(datetime.datetime.now()) + " MECAT correction done!"

def makePrimerForLimaAndIsoseq(primer, isoseqPrimer, limaPrimer):
    records = list(SeqIO.parse(primer, "fasta"))
    if len(records) < 2:
        sys.exit("The primer is not correct!")
    isoseqOut = open(isoseqPrimer, "w")
    limaOut = open(limaPrimer, "w")
    for i in range(len(records)):
        if i == 0:
            print >>isoseqOut, ">F0"
            print >>isoseqOut, records[i].seq
            print >>limaOut, ">5p"
            print >>limaOut, records[i].seq
        else:
            print >>isoseqOut, ">R0"
            print >>isoseqOut, records[i].seq
            print >>limaOut, ">3p"
            print >>limaOut, records[i].seq
    isoseqOut.close()
    limaOut.close()

def smrtOldMethods(ccsBam, paramObj=None, ccsParams=None):
    tgs_dir, ngs_dir, pac_plat, projectName = paramObj.tgs_dir, paramObj.ngs_dir, paramObj.pac_plat, paramObj.project
    primer = os.path.join(tgs_dir, "primer.fasta")
    # 1.2 Extract raw reads
    # print str(datetime.datetime.now()) + " Extract reads..."
    # inputBam = os.path.join(tgs_dir, "input.subreads.bam")
    # bamProcess(inputBam, "input.rawsubreads.fastq", bam2Fastq=True)
    # cmd = "samtools fastq --threads {} {} > input.rawsubreads.fastq 2> " \
    #       "log/samtools.rawsubreads.fastq.log".format(thread, inputBam)
    # cmd = "fqNameMapping.pl -p {} input.rawsubreads.fastq >input.rawsubreads.fq 2>rawsubreads.idmapping & ".format(
    #     projectName)
    # print str(datetime.datetime.now()) + " Extract rawreads and subreads done!"

    # 1.3 Get CCS.fa
    print str(datetime.datetime.now()) + " Extract CCS..."
    # cmd = "ccs {} CCS.bam -j {} --noPolish --minPasses {} --minLength {} --minPredictedAccuracy {} --minReadScore " \
    #       "{}".format(inputBam, thread, ccsParams.min_pass, ccsParams.min_subread_length,
    #                   ccsParams.min_predicted_accuracy, ccsParams.min_read_score)
    bamProcess(ccsBam, outPrefix="CCS", name2pass=True)
    # cmd = "samtools view CCS.bam | cut -f1,13 | sed 's/np:i://' >name2pass.tsv"
    # cmd = "samtools fastq --threads {} CCS.bam 2>log/samtools.ccs.fastq.log | fqNameMapping.pl -p {}/ >CCS.fq 2>CCS.idmapping".format(
    #     thread, projectName)
    # cmd = '''fastqToFa CCS.fq CCS.fa'''
    print str(datetime.datetime.now()) + " Extract CCS done!"

    # 2. Get full-length(FL) reads and pre-mapping evaluation
    # 2.1 Get FL reads
    primerSearchWin, primerSearchMinScore = 150, 8
    print str(datetime.datetime.now()) + " Get full-length (FL) reads..."
    hmmerWrapperCommonParams = "--primer_search_window {} --min-score {} --cpus " \
                               "{} --change-seqid --input_filename CCS.fa --output_filename FL.fa --directory " \
                               "hmmer/ --must-see-polyA".format(primerSearchWin, primerSearchMinScore, thread)
    cmd = "(time hmmer_wrapper.py {} --primer_filename {} 2>log/hmmer_wrapper.log) 2>log/hmmer_wrapper.time".format(
        hmmerWrapperCommonParams, primer)
    cmd = "summarize_primer_info.py FL.fa.primer_info.txt >log/summarize_primer_info.log"
    os.rename("FL.fa.primer_info.txt", "primer_info.tsv")
    # cmd = "mv FL.fa.primer_info.txt primer_info.tsv"

    # 2.2 Identify chimera
    cmd = "(time chimera_finder.py --min_dist-from_end 50 --cpus {} --primer_filename {} --input_filename " \
          "FL.fa --directory hmmer_chimera/ 2>log/chimera_finder.FL.log) 2>log/chimera_finder.time".format(thread,
                                                                                                           primer)
    os.rename("FL.fa.non_chimera.fa", "FLnoC.fa")
    os.rename("FL.fa.is_chimera.fa", "FLC.fa")
    # cmd = "mv FL.fa.non_chimera.fa FLnoC.fa & mv FL.fa.is_chimera.fa FLC.fa"

    # 2.3 Trim polyA
    windowSize = 50
    fraction = 0.65
    cmd = "faTrimPA.pl -w {} -f {} FLnoC.fa 2>log/faTrimPA.log | " \
          "faFilter -minSize={} stdin paTrimmed.fa".format(windowSize, fraction, ccsParams.min_ccs_length)
    getPAtailLength("log/faTrimPA.log", "primer_info.tsv", "paTailLength.tsv")
    # cmd = '''( awk '$3<$2' log/faTrimPA.log | sed 's/^>//' |
    #                 skyjoin - primer_info.tsv 2>/dev/null | awk '$10>$9{print $1"\t"$10-$9+$2-$3}'
    #            awk '$3<$2' log/faTrimPA.log |
    #                 sed 's/^>//' | filter.pl -o /dev/stdin <(grep -v '^ID' primer_info.tsv) |
    #                 awk '$7!="NA"&&$8!="NA"&&$8-$7>0{print $1"\t"$8-$7}') | sort -n >paTailLength.tsv
    #             '''
    cmd = "cut -f2 paTailLength.tsv | hist.R -x='Tail Length' -y=Density -b=1 -d -x1=0 -x2=100 -p=paTailLength.pdf 2>/dev/null"
    print str(datetime.datetime.now()) + " Get full-length (FL) reads done!"

    # 2.4 Pre-mapping evaluation
    print str(datetime.datetime.now()) + " Pre-mapping evaluation..."
    os.makedirs("evaluation")
    os.chdir("evaluation")
    cmd = "faCount ../paTrimmed.fa >paTrimmed.faCount"
    gc_out = open("GC_of_reads.log", "w")
    with open("paTrimmed.faCount") as faCount:
        for line in faCount:
            if line.startswith("#") or "total" in line: continue
            faCountList = line.strip().split("\t")
            print >> gc_out, "{}\t{}".format(faCountList[0],
                                             (int(faCountList[3]) + int(faCountList[4])) * 100 / float(faCountList[1]))
    gc_out.close()
    # cmd = '''grep -v '^#' paTrimmed.faCount | grep -v total | awk '{print $1"\t"($4+$5)/$2*100}' >GC_of_reads.log'''
    cmd = '''cut -f2 GC_of_reads.log | distrCurve.R -d -m='GC Content of Reads' -x='Binned GC%' -y='Fraction of Reads' -v=50 -p=GC_of_reads.pdf 2>/dev/null'''
    cmd = "GC_across_read.pl -i 20 ../paTrimmed.fa >GC_across_read.log"
    cmd = '''cut -f2- GC_across_read.log | box.R -stack -nJ -ho=50 -m='GC Content across Reads' -x=Interval -y=GC% -oS=0.5 -w=11 -p=GC_across_read.pdf 2>/dev/null'''
    cmd = '''BQ_across_read.pl ../CCS.fq | tee BQ_across_read.log | cut -f2- | box.R -nJ -stack -m='Quality across Reads' -x=Interval -y=Quality -w=11 -p=BQ_across_read.pdf 2>/dev/null'''
    print str(datetime.datetime.now()) + " Pre-mapping evaluation done!"

def isoseq3Corr(inputBam, subreadsXML, paramObj=None, refParams=None, hybrid=False):
    print str(datetime.datetime.now()) + " IsoSeq3 correction..."
    prevDir = os.getcwd()
    if os.path.exists("isoseq3/primer"):
        shutil.rmtree("isoseq3/primer")
    os.makedirs("isoseq3/primer")
    os.chdir("isoseq3")
    shutil.copy("{}/primer.fasta".format(paramObj.tgs_dir), "primer/")
    # os.symlink("{}/primer.fasta".format(paramObj.tgs_dir), "primer/")
    # cmd = "cp {}/primer.fasta primer/".format(paramObj.tgs_dir)
    isoseqPrimer = "primerIsoseq.fa"
    limaPrimer = "primerLima.fa"
    makePrimerForLimaAndIsoseq("primer/primer.fasta", isoseqPrimer, limaPrimer)
    # cmd = "makePrimerForLimaAndIsoseq.py primer/primer.fasta 1>primerIsoseq.fa 2>primerLima.fa"
    cmd = "lima --isoseq --dump-clips --no-pbi -j {} ../{} primerLima.fa demux.bam".format(thread, inputBam)
    # print cmd
    subprocess.call(cmd, shell=True)
    cmd = "isoseq3 refine --require-polya demux.5p--3p.bam primerLima.fa flnc.bam"
    # print cmd
    subprocess.call(cmd, shell=True)
    cmd = "isoseq3 cluster -j {} flnc.bam unpolish.flnc.bam".format(thread)
    # print cmd
    subprocess.call(cmd, shell=True)
    cmd = "isoseq3 polish -j {} unpolish.flnc.bam ../{} polished.flnc.bam".format(thread, subreadsXML)
    # print cmd
    subprocess.call(cmd, shell=True)
    cmd = "samtools fastq -@ {} polished.flnc.bam > polished.flnc.fq".format(thread)
    # print cmd
    subprocess.call(cmd, shell=True)
    if hybrid:
        os.chdir(prevDir)
        print str(datetime.datetime.now()) + " IsoSeq3 correction done!"
        return os.path.join(prevDir, "isoseq3", "polished.flnc.fq")
    else:
        cmd = "minimap2 -ax splice:hq -uf {} -t {} --MD polished.flnc.fq > aln.sam".format(refParams.ref_mm2_index, thread)
        subprocess.call(cmd, shell=True)
        os.chdir(prevDir)
        os.symlink(os.path.join(prevDir, "isoseq3", "aln.sam"), "aln.sam")
        print str(datetime.datetime.now()) + " IsoSeq3 correction done!"

def correctWithDifferentTools(seqNeedCorrect, paramObj=None, refParams=None, mecatParams=None, inputBam=None,
                              subreads=None, ngs_dir=None, tgs_plat="pacbio"):
    if not os.path.exists(refParams.ref_mm2_index):
        cmd = "minimap2 -d {} {}".format(refParams.ref_mm2_index, refParams.ref_genome)
        subprocess.call(cmd, shell=True)
        # print cmd
    if tgs_plat == "pacbio" and ngs_dir:
        isoseqCorrected = isoseq3Corr(inputBam, subreads, paramObj=paramObj, refParams=refParams, hybrid=True)
        proovreadCorr(isoseqCorrected, paramObj=paramObj, refParams=refParams, tgs_plat="pacbio")
    elif tgs_plat == "pacbio" and ngs_dir == None:
        isoseq3Corr(inputBam, subreads, paramObj=paramObj, refParams=refParams)
    elif tgs_plat == "nanopore" and ngs_dir:
        mecatCorrected = mecatCorr(seqNeedCorrect, refParams=refParams, mecatParams=mecatParams, hybrid=True)
        proovreadCorr(mecatCorrected, paramObj=paramObj, refParams=refParams, tgs_plat="nanopore")
    elif tgs_plat == "nanopore" and ngs_dir == None:
        mecatCorr(seqNeedCorrect, refParams=refParams, mecatParams=mecatParams)
    elif ngs_dir:
        print str(datetime.datetime.now()) + " No TGS reads provided, only NGS analysis will be conducted!"

    cmd = "samtools view -bS aln.sam | samtools sort -@ 8 -o aln.sorted.bam && samtools index aln.sorted.bam"
    subprocess.call(cmd, shell=True)

def collapse(tgs_plat="pacbio", refParams=None):
    print str(datetime.datetime.now()) + " Collapse with cDNA_cupcake..."
    prevDir = os.getcwd()
    if os.path.exists("collapse"):
        shutil.rmtree("collapse")
    os.makedirs("collapse")
    os.chdir("collapse")
    os.symlink("../sqanti/sqanti_artifact_removed.fa", "sqanti_artifact_removed.fa")
    os.symlink("../mapping/processed.bed12+", "processed.bed12+")
    if tgs_plat == "pacbio":
        cmd = "minimap2 -ax splice:hq -uf --secondary=no -t {} {} sqanti_artifact_removed.fa | samtool view -bS | samtool sort | samtool view -h -O sam > artifact_removed.sorted.sam".format(thread, refParams.ref_mm2_index)
    else:
        cmd = "minimap2 -ax splice -uf --secondary=no -t {} -k14 {} sqanti_artifact_removed.fa | samtool view -bS | samtool sort | samtool view -h -O sam > artifact_removed.sorted.sam".format(thread, refParams.ref_mm2_index)
    subprocess.call(cmd, shell=True)

    cmd = "collapse_isoforms_by_sam.py --input sqanti_artifact_removed.fa -s artifact_removed.sorted.sam --dun-merge-5-shorter -o tofu"
    subprocess.call(cmd, shell=True)
    cmd = '''echo -e "pbid\tcount_fl" > tofu.collapsed.abundance.txt'''
    subprocess.call(cmd, shell=True)
    cmd = '''paste <(awk '{print $1}' tofu.collapsed.group.txt) <(awk -F , '{print NF}' tofu.collapsed.group.txt ) >> tofu.collapsed.abundance.txt'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "grep -v '@' artifact_removed.sorted.sam | cut -f1 > artifact_removed.read.lst"
    subprocess.call(cmd, shell=True)
    cmd = "filter.pl -o tofu.ignored_ids.txt artifact_removed.read.lst | tee artifact_removed.ignored_id_removed.read.lst | seqkit grep sqanti_artifact_removed.fa -f - -w 0 > sqanti_artifact_removed.ignore_id_removed.fa"
    subprocess.call(cmd, shell=True)
    cmd = "filter.pl -o artifact_removed.ignored_id_removed.read.lst processed.bed12+ -2 4 -m i > processed.ignore_id_removed.bed12+"
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Collapse with cDNA_cupcake done!"

def buildSpliceGrapherModel(refParams=None):
    print str(datetime.datetime.now()) + " Start SpliceGrapher build classification modeling..."
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
    print str(datetime.datetime.now()) + " SpliceGrapher build classification modeling done!"

class Bed6(object):
    def __init__(self, line=None):
        if line:
            self.record = line.strip("\n").split("\t")
            self.chrom, self.start, self.end, self.name, self.score, self.strand = self.record[:6]
            self.start = int(self.start)
            self.end = int(self.end)
            self.score = int(self.score) if self.score != '.' else self.score
        else:
            self.empty()

    def empty(self):
        (self.chrom, self.start, self.end, self.name, self.score, self.strand) = ('', 0, 0, '', 0, '')

    def __str__(self):
        return "%s:%d-%d:%s,%s" % (self.chrom, self.start, self.end, self.strand, str(self.score))

class Bed12(object):
    def __init__(self, line=None):
        if line:
            self.record = line.strip("\n").split("\t")
            (self.chrom, self.start, self.end, self.name,
             self.score, self.strand, self.thickStart, self.thickEnd,
             self.rgb, self.blockCount, self.blockSizes,
             self.blockStarts) = self.record[:12]
            self.start = int(self.start)
            self.end = int(self.end)
            self.score = int(self.score)
            self.thickStart = int(self.thickStart)
            self.thickEnd = int(self.thickEnd)
            self.blockCount = int(self.blockCount)

            self.blockSizes = map(int, self.blockSizes.strip(",").split(","))
            self.blockStarts = map(int, self.blockStarts.strip(",").split(","))

            self.exonStarts = map(lambda x: x+self.start, self.blockStarts)
            self.exonEnds = map(lambda x: self.exonStarts[x]+self.blockSizes[x], range(self.blockCount))
            self.exonList = map(lambda x: [self.exonStarts[x], self.exonEnds[x]], range(self.blockCount))
            self.minpos = min(self.exonStarts)
            self.maxpos = max(self.exonEnds)
        else:
            self.empty()

    def empty(self):
        self.chrom = ""
        self.start, self.end = 0, 0
        self.name = ""
        self.score = 0
        self.strand = ""
        self.thickStart = self.thickEnd = 0
        self.rgb = "0,0,0"
        self.blockCount = 0
        self.blockSizes, self.blockStarts = [], []
        self.exonStarts, self.exonEnds = [], []
        self.exonList = []
        self.minpos, self.maxpos = 0, 0

    def __len__(self):
        return sum(self.blockSizes)

    def __str__(self):
        fields = [self.chrom, str(self.start), str(self.end),
                  self.name, str(self.score), self.strand, str(self.thickStart),
                  str(self.thickEnd), self.rgb, str(self.blockCount)]
        blockSizesLine = ','.join(repr(i) for i in self.blockSizes)
        blockStartsLine = ','.join(repr(i) for i in self.blockStarts)
        fields.extend([blockSizesLine, blockStartsLine])
        return "\t".join(fields)

    __repr__ = __str__

class Bed6Plus(Bed6):
    def __init__(self, line=None):
        if line:
            Bed6.__init__(self, line)
            self.otherList = self.record[6:]

    def __str__(self):
        return Bed6.__str__(self) + "\t" + "\t".join(self.otherList)

class Bed12Plus(Bed12):
    def __init__(self, line=None):
        if line:
            Bed12.__init__(self, line)
            self.otherList = self.record[12:]

    def __str__(self):
        return Bed12.__str__(self) + "\t" + "\t".join(self.otherList)

class BedFile(object):
    def __init__(self, bedFile, type=None):
        self.bedFile = bedFile
        self.reads = self.getReadsInfo(type)

    def getReadsInfo(self, type):
        readsDict = {}
        with open(self.bedFile) as f:
            for line in f:
                if type == "bed12":
                    b = Bed12(line)
                    readsDict.__setitem__(b.name, b)
                elif type == "bed12+":
                    b = Bed12Plus(line)
                    readsDict.__setitem__(b.name, b)
                elif type == "bed6":
                    b = Bed6(line)
                    readsDict.__setitem__(b.name, b)
                elif type == "bed6+":
                    b = Bed6Plus(line)
                    readsDict.__setitem__(b.name, b)
        return readsDict

def filterAndRefineAlign(refParams=None):
    print str(datetime.datetime.now()) + " Filter and refine alignments..."
    prevDir = os.getcwd()
    if os.path.exists("mapping"):
        shutil.rmtree("mapping")
    os.makedirs("mapping")
    os.chdir("mapping")
    buildSpliceGrapherModel(refParams=refParams)
    # print str(datetime.datetime.now()) + " Start SpliceGrapher build classification modeling..."
    # cmd = "generate_splice_site_data.py -f {} -m {} -r spliceSiteFreq.txt 1>/dev/null".format(refParams.ref_genome, refParams.ref_gtf)
    # juncFreqDict = {}
    # with open("spliceSiteFreq.txt") as f:
    #     lineList = f.readlines()
    #     juncIndex = lineList.index("Breakdown of splice junctions:")
    #     sumFreq = 0
    #     for i in range(juncIndex+1, len(lineList)):
    #         junc, info = re.split(":", lineList[i].strip())
    #         donor, acceptor = junc.strip().split("-")
    #         freq = float(re.search("\((.*)%\)", info).groups()[0])
    #         sumFreq += freq
    #         juncFreqDict.__setitem__("donor", donor)
    #         juncFreqDict.__setitem__("acceptor", acceptor)
    #         if sumFreq >= 98: break
    #
    # donorClass = ",".join(juncFreqDict["donor"])
    # acceptorClass = ",".join(juncFreqDict["acceptor"])
    # cmd = "build_classifiers.py -m {} -f {} -n 2000 -d {} -a {}".format(refParams.ref_gtf, refParams.ref_genome,
    #                                                                     donorClass, acceptorClass)
    # print str(datetime.datetime.now()) + " SpliceGrapher build classification modeling done!"
    cmd = '''samAddTag.pl --checkHardClip --coverage --identity --unmapped unmapped.sam ../aln.sam 2>lengthInconsistent.sam | samtools view -buS - | samtools sort -m 4G - >mapped.sorted.bam'''
    subprocess.call(cmd, shell=True)
    cmd = "sam2bed.pl -t CV,ID mapped.sorted.bam >mapped.bed12+"
    subprocess.call(cmd, shell=True)
    cmd = "readsFilter.pl -r 0.8 mapped.bed12+ 2>discarded.bed12+ >uniq.bed12+"
    subprocess.call(cmd, shell=True)
    cmd = r'''(samtools view -H mapped.sorted.bam;samtools view mapped.sorted.bam | filter.pl -o <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2+1,$4,$5}' uniq.bed12+) -1 1,2,3,4 -2 3,4,1,5 -m i) > uniq.sam'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    print str(datetime.datetime.now()) + " Start filtering sam file with model built by SpliceGrapher..."
    cmd = "sam_filter.py uniq.sam classifiers.zip -f {} -o filtered.uniq.sam".format(refParams.ref_genome)
    subprocess.call(cmd, shell=True)
    # pysam.sort("-m", "4G", "-@", "5", "-o", "uniq1.sorted.bam", "filtered.uniq.sorted.sam")
    cmd = "samtools sort -m 4G filtered.uniq.sam >uniq.sorted.bam"
    subprocess.call(cmd, shell=True)
    print str(datetime.datetime.now()) + " Filtering sam file with model build by SpliceGrapher done!"

    cmd = "strandAdjust.pl -f {} -g {} -c 0.8 -s 2 uniq.bed12+ >strandAdjusted.bed12+".format(refParams.ref_genome, refParams.ref_gpe)
    subprocess.call(cmd, shell=True)

    # confirm strand
    uniqBed = BedFile("uniq.bed12+", type="bed12+")
    strandAdjustedBed = BedFile("strandAdjusted.bed12+", type="bed12+")
    strandConfirm = open("strandConfirm.bed12+", "w")

    for read in uniqBed.reads:
        if uniqBed.reads[read].strand != strandAdjustedBed.reads[read].strand: continue
        print >>strandConfirm, uniqBed.reads[read]
    strandConfirm.close()
    # cmd = "paste <(cut -f4,6 uniq.bed12+) <(cut -f6 strandAdjusted.bed12+) | awk '$2!=$3{print $1}' >strandConflict.readName"
    # cmd = "filter.pl -o strandConflict.readName -2 4 uniq.bed12+ >strandConfirm.bed12+"

    juncScoringParams = "-r {} exonRealign.bed12+".format(refParams.ref_gpe)
    # if paramObj.ngs_dir:
    #     cmd = "removeMisAlignExon.pl -g ../../raw.gpe -j $junction strandConfirm.bed12+ | fillMissExon.pl -g ../../raw.gpe -b -j $junction >exonRealign.bed12+"
    #     cmd = '''juncScoringParams="{} -j $junction"'''.format(juncScoringParams)
    # else:
    cmd = "removeMisAlignExon.pl -g {} -j /dev/null strandConfirm.bed12+ | fillMissExon.pl -g {} -b -j /dev/null >exonRealign.bed12+".format(refParams.ref_gpe, refParams.ref_gpe)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "juncConsensus.pl -s <(juncScoring.pl {}) exonRealign.bed12+ >consensus.bed12+".format(juncScoringParams)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    # cmd = "bed2gpe.pl consensus.bed12+ | gpeFeature.pl --exon | filter.pl -o /dev/stdin -1 4 -2 4 consensus.bed12+ >noGap.bed12+"
    # subprocess.call(cmd, shell=True)
    # cmd = '''dnaOrInternalPrimingContFilter.pl -b {} -t <(sed 's/^@//' ../paTailLength.tsv) -r dnaOrInternalPrimingContFilter.log consensus.bed12+ >deCont.bed12+ 2>dnaCont.bed12+'''.format(refParams.ref2bit)
    # subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = r'''(samtools view -F 0x10 uniq.sorted.bam | perl -ne 'print if (split "\t")[5] =~ /(\d+)S$/ && $1 >30';samtools view -f 0x10 uniq.sorted.bam | perl -ne 'print if (split "\t")[5] =~ /^(\d+)S/ && $1 >30') | filter.pl -o /dev/stdin -2 4 consensus.bed12+ >processed.bed12+'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "samtools fasta ../aln.sam 2>/dev/null | seqkit grep - -f <(cut -f4 processed.bed12+) -w 0 >processed.fa"
    subprocess.call(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE)

    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Filter and refine results done!"

def removeArtifacts(refParams=None):
    print str(datetime.datetime.now()) + " Remove artifact isoforms..."
    prevDir = os.getcwd()
    if os.path.exists("sqanti"):
        shutil.rmtree("sqanti")
    os.makedirs("sqanti")
    os.chdir("sqanti")
    os.symlink("../mapping/processed.bed12+", "processed.bed12+")
    os.symlink("../mapping/processed.fa", "processed.fa")
    cmd = "bedToGenePred processed.bed12+ processed.gp"
    subprocess.call(cmd, shell=True)
    cmd = "genePredToGtf file processed.gp processed.gtf"
    subprocess.call(cmd, shell=True)
    cmd = "sqanti_qc.py processed.gtf {} {} -g -d ./ -o sqanti".format(refParams.ref_gtf, refParams.ref_genome)
    subprocess.call(cmd, shell=True)
    cmd = "sqanti_filter.py sqanti_classification.txt"
    subprocess.call(cmd, shell=True)
    cmd = "seqkit grep processed.fa -f Filter_out/sqanti_classification.txt_curatedTranscriptome.txt -w 0 >sqanti_artifact_removed.fa"
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Remove artifact isoforms done!"

def getSRfileNames(ngs_dir):
    fileDict = {}
    for i in os.listdir(ngs_dir):
        if os.path.isfile(i):
            if len(re.compile(".r1.", re.IGNORECASE).split(i)) > 1:
                filePrex = re.compile(".r1.", re.IGNORECASE).split(i)[0]
                fileDict[filePrex].__setitem__("r1", i)
            if len(re.compile(".r2.", re.IGNORECASE).split(i)) > 1:
                filePrex = re.compile(".r2.", re.IGNORECASE).split(i)[0]
                fileDict[filePrex].__setitem__("r2", i)
    return fileDict

def bamProcess(bamFile, outPrefix="out", name2pass=False, bam2Fastx=False):
    bam = pysam.AlignmentFile(bamFile, "rb", check_sq=False)
    outFq = open(outPrefix + ".fq", "w")
    outFa = open(outPrefix + ".fa", "w")
    name2passOut = open(outPrefix + ".name2pass.tsv", "w")
    for line in bam:
        if bam2Fastx:
            print >>outFq, "@{}".format(line.query_name)
            print >>outFq, "{}".format(line.query_sequence)
            print >>outFq, "+"
            print >>outFq, "{}".format(line.query_qualities)
            print >>outFa, ">{}".format(line.query_name)
            print >>outFa, "{}".format(line.query_sequence)
        if name2pass:
            print >>name2passOut, "{}\t{}".format(line.query_name, line.get_tag("np"))

    outFa.close()
    outFq.close()
    name2passOut.close()

def getPAtailLength(faTrimPA, primerInfo, paTailOut):
    trimPADict = {}
    primerDict = {}
    paTailLength = {}
    out = open(paTailOut, "w")
    with open(faTrimPA) as f1:
        for line1 in f1:
            infoList = line1.strip().split("\t")
            read = infoList[0].split(">")[1]
            trimPADict[read] = map(int, infoList[1:])
    with open(primerInfo) as f2:
        for line2 in f2:
            if line2.startswith("ID"): continue
            infoList = line2.strip().split("\t")
            primerDict[infoList[0]] = infoList[2:]
    for read in primerDict:
        if read in trimPADict and trimPADict[read][1] < trimPADict[read][0]:
            paTailLength[read] = trimPADict[read][0] - trimPADict[read][1] + int(primerDict[read][5]) - int(primerDict[read][4])
            print >>out, "{}\t{}".format(read, paTailLength[read])
        else:
            if primerDict[read][5] != "NA" and primerDict[read][4] != "NA" and int(primerDict[read][5]) > int(primerDict[read][4]):
                paTailLength[read] = int(primerDict[read][5]) - int(primerDict[read][4])
                print >>out, "{}\t{}".format(read, paTailLength[read])
    out.close()

def processReadsFromPacbio(paramObj=None, refParams=None, ccsParams=None):
    tgs_dir, ngs_dir, pac_plat, projectName = paramObj.tgs_dir, paramObj.ngs_dir, paramObj.pac_plat, paramObj.project
    prevDir = os.getcwd()
    # 1. Extract reads
    # 1.1 Determine the sequencing strategy
    os.mkdir("log")
    primer = os.path.join(tgs_dir, "primer.fasta")
    print str(datetime.datetime.now()) + " Extract reads..."
    if pac_plat.lower() == "rsii":
        print str(datetime.datetime.now()) + " The reads are from RSII"
        cmd = "find {} -name '*.bax.h5' > input.fofn".format(tgs_dir)
        subprocess.call(cmd, shell=True)
        cmd = "bax2bam -f input.fofn -o input"
        subprocess.call(cmd, shell=True)

        inputBam = subreads = "input.subreads.bam"
        print str(datetime.datetime.now()) + " Extract CCS reads from bax.h5 files..."
        cmd = "ccs {} CCS.bam -j {} --noPolish --minPasses {} --minLength {} --minPredictedAccuracy {} --minReadScore " \
              "{}".format(inputBam, thread, ccsParams.min_pass, ccsParams.min_subread_length,
                          ccsParams.min_predicted_accuracy, ccsParams.min_read_score)
        subprocess.call(cmd, shell=True)
        print str(datetime.datetime.now()) + " Extract CCS reads from bax.h5 files done!"
        # primer = "primer.fasta"
    elif pac_plat.lower() == "sequel":
        print str(datetime.datetime.now()) + " The reads are from Sequel"
        xmlFiles = os.popen("find -L {} -name '*.subreadset.xml'".format(tgs_dir))
        xmlFileStr = " ".join(xmlFiles)
        # os.makedirs("{}/{}".format(refParams.out_dir, projectName))
        subreads = "input.subreadset.xml"
        cmd = "dataset create --type SubreadSet --force --name {} {} {}".format(projectName, subreads, xmlFileStr)
        subprocess.call(cmd, shell=True)
        print str(datetime.datetime.now()) + " Extract CCS reads from subreads files..."
        cmd = "ccs {} CCS.bam -j {} --noPolish --minPasses {} --minLength {} --minPredictedAccuracy {} --minReadScore " \
              "{}".format(subreads, thread, ccsParams.min_pass, ccsParams.min_subread_length,
                          ccsParams.min_predicted_accuracy, ccsParams.min_read_score)
        subprocess.call(cmd, shell=True)
        print str(datetime.datetime.now()) + " Extract CCS reads from subreads files done!"
        # subprocess.call(
        #     "ccs {} CCS.bam -j {} --noPolish --minPasses {} --minLength {} --minPredictedAccuracy {} ".format(
        #         outXML, thread, ccsParams.min_pass, ccsParams.min_subread_length, ccsParams.min_predicted_accuracy), shell=True)
        # return
    else:
        raise Exception("You should input the right PacBio sequencing strategy!!")

    # 1.2 Extract raw reads
    # print str(datetime.datetime.now()) + " Extract reads..."
    # inputBam = os.path.join(tgs_dir, "input.subreads.bam")
    # bamProcess(inputBam, "input.rawsubreads.fastq", bam2Fastq=True)
    # cmd = "samtools fastq --threads {} {} > input.rawsubreads.fastq 2> " \
    #       "log/samtools.rawsubreads.fastq.log".format(thread, inputBam)
    # cmd = "fqNameMapping.pl -p {} input.rawsubreads.fastq >input.rawsubreads.fq 2>rawsubreads.idmapping & ".format(
    #     projectName)
    # print str(datetime.datetime.now()) + " Extract rawreads and subreads done!"

    # 1.3 Get CCS.fa
    print str(datetime.datetime.now()) + " Extract CCS..."
    # cmd = "ccs {} CCS.bam -j {} --noPolish --minPasses {} --minLength {} --minPredictedAccuracy {} --minReadScore " \
    #       "{}".format(inputBam, thread, ccsParams.min_pass, ccsParams.min_subread_length,
    #                   ccsParams.min_predicted_accuracy, ccsParams.min_read_score)
    bamProcess("CCS.bam", outPrefix="CCS", name2pass=True, bam2Fastx=True)
    # cmd = "samtools view CCS.bam | cut -f1,13 | sed 's/np:i://' >name2pass.tsv"
    # cmd = "samtools fastq --threads {} CCS.bam 2>log/samtools.ccs.fastq.log | fqNameMapping.pl -p {}/ >CCS.fq 2>CCS.idmapping".format(
    #     thread, projectName)
    # cmd = '''fastqToFa CCS.fq CCS.fa'''
    print str(datetime.datetime.now()) + " Extract CCS done!"

    # 2. Get full-length(FL) reads and pre-mapping evaluation
    # 2.1 Get FL reads
    primerSearchWin, primerSearchMinScore = 150, 8
    print str(datetime.datetime.now()) + " Get full-length (FL) reads..."
    hmmerWrapperCommonParams = "--primer_search_window {} --min-score {} --cpus " \
                               "{} --change-seqid --input_filename CCS.fa --output_filename FL.fa --directory " \
                               "hmmer/ --must-see-polyA".format(primerSearchWin, primerSearchMinScore, thread)
    cmd = "(time hmmer_wrapper.py {} --primer_filename {} 2>log/hmmer_wrapper.log) 2>log/hmmer_wrapper.time".format(
        hmmerWrapperCommonParams, primer)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    # subprocess.call("(time hmmer_wrapper.py {} --primer_filename {} 2>log/hmmer_wrapper.log) 2>log/hmmer_wrapper.time".format(
    #     hmmerWrapperCommonParams, primer), shell=True)
    subprocess.call("summarize_primer_info.py FL.fa.primer_info.txt >log/summarize_primer_info.log", shell=True)
    # cmd = "summarize_primer_info.py FL.fa.primer_info.txt >log/summarize_primer_info.log"
    os.rename("FL.fa.primer_info.txt", "primer_info.tsv")
    # cmd = "mv FL.fa.primer_info.txt primer_info.tsv"

    # 2.2 Identify chimera
    cmd = "(time chimera_finder.py --min_dist-from_end 50 --cpus {} --primer_filename {} --input_filename FL.fa " \
          "--directory hmmer_chimera/ 2>log/chimera_finder.FL.log) 2>log/chimera_finder.time".format(
        thread, primer)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    os.rename("FL.fa.non_chimera.fa", "FLnoC.fa")
    os.rename("FL.fa.is_chimera.fa", "FLC.fa")

    # 2.3 Trim polyA
    windowSize = 50
    fraction = 0.65
    cmd = "faTrimPA.pl -w {} -f {} FLnoC.fa 2>log/faTrimPA.log | " \
          "faFilter -minSize={} stdin paTrimmed.fa".format(windowSize, fraction, ccsParams.min_ccs_length)
    subprocess.call(cmd, shell=True)
    getPAtailLength("log/faTrimPA.log", "primer_info.tsv", "paTailLength.tsv")
    # cmd = '''( awk '$3<$2' log/faTrimPA.log | sed 's/^>//' |
    #                 skyjoin - primer_info.tsv 2>/dev/null | awk '$10>$9{print $1"\t"$10-$9+$2-$3}'
    #            awk '$3<$2' log/faTrimPA.log |
    #                 sed 's/^>//' | filter.pl -o /dev/stdin <(grep -v '^ID' primer_info.tsv) |
    #                 awk '$7!="NA"&&$8!="NA"&&$8-$7>0{print $1"\t"$8-$7}') | sort -n >paTailLength.tsv
    #             '''
    cmd = "cut -f2 paTailLength.tsv | hist.R -x='Tail Length' -y=Density -b=1 -d -x1=0 -x2=100 -p=paTailLength.pdf 2>/dev/null"
    subprocess.call(cmd, shell=True)
    print str(datetime.datetime.now()) + " Get full-length (FL) reads done!"

    # 2.4 Pre-mapping evaluation
    print str(datetime.datetime.now()) + " Pre-mapping evaluation..."
    os.makedirs("evaluation")
    os.chdir("evaluation")
    cmd = "faCount ../paTrimmed.fa >paTrimmed.faCount"
    subprocess.call(cmd, shell=True)
    gc_out = open("GC_of_reads.log", "w")
    with open("paTrimmed.faCount") as faCount:
        for line in faCount:
            if line.startswith("#") or "total" in line: continue
            faCountList = line.strip().split("\t")
            print >>gc_out, "{}\t{}".format(faCountList[0], (int(faCountList[3])+int(faCountList[4]))*100/float(faCountList[1]))
    gc_out.close()
    # cmd = '''grep -v '^#' paTrimmed.faCount | grep -v total | awk '{print $1"\t"($4+$5)/$2*100}' >GC_of_reads.log'''
    cmd = '''cut -f2 GC_of_reads.log | distrCurve.R -d -m='GC Content of Reads' -x='Binned GC%' -y='Fraction of Reads' -v=50 -p=GC_of_reads.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    cmd = "GC_across_read.pl -i 20 ../paTrimmed.fa >GC_across_read.log"
    subprocess.call(cmd, shell=True)
    cmd = '''cut -f2- GC_across_read.log | box.R -stack -nJ -ho=50 -m='GC Content across Reads' -x=Interval -y=GC% -oS=0.5 -w=11 -p=GC_across_read.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    cmd = '''BQ_across_read.pl ../CCS.fq | tee BQ_across_read.log | cut -f2- | box.R -nJ -stack -m='Quality across Reads' -x=Interval -y=Quality -w=11 -p=BQ_across_read.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Pre-mapping evaluation done!"

    # 3 Process with different correction tools
    seqNeedCorrect = "paTrimmed.fa"
    inputBam = "CCS.bam"
    return inputBam, seqNeedCorrect, subreads
    # correctWithDifferentTools(seqNeedCorrect, paramObj=paramObj, refParams=refParams, inputBam=inputBam,
    #                           ngs_dir=ngs_dir, tgs_plat="pacbio")

def preMappingEvaluation():
    print str(datetime.datetime.now()) + " Pre-mapping evaluation..."
    os.makedirs("evaluation")
    os.chdir("evaluation")
    cmd = "faCount ../paTrimmed.fa >paTrimmed.faCount"
    cmd = '''grep -v '^#' paTrimmed.faCount | grep -v total | awk '{print $1"\t"($4+$5)/$2*100}' >GC_of_reads.log'''
    cmd = '''cut -f2 GC_of_reads.log | distrCurve.R -d -m='GC Content of Reads' -x='Binned GC%' -y='Fraction of Reads' -v=50 -p=GC_of_reads.pdf 2>/dev/null'''
    cmd = "GC_across_read.pl -i 20 ../paTrimmed.fa >GC_across_read.log"
    cmd = '''cut -f2- GC_across_read.log | box.R -stack -nJ -ho=50 -m='GC Content across Reads' -x=Interval -y=GC% -oS=0.5 -w=11 -p=GC_across_read.pdf 2>/dev/null'''
    cmd = '''BQ_across_read.pl ../CCS.fq | tee BQ_across_read.log | cut -f2- | box.R -nJ -stack -m='Quality across Reads' -x=Interval -y=Quality -w=11 -p=BQ_across_read.pdf 2>/dev/null'''
    print str(datetime.datetime.now()) + " Pre-mapping evaluation done!"

def processReadsFromNanopore():
    pass

def degOnlySR(bamList, paramObj=None):
    if paramObj.ngs_sample_cond == None:
        sampleCond = open("sampleCond.txt", "w")
        print >>sampleCond, "sample_id\tconditions"
        for i in bamList:
            print >>sampleCond, "{}\t{}".format(i, i)
        sampleCond.close()
        paramObj.ngs_sample_cond = os.path.join(os.getcwd(), "sampleCond.txt")
    cmd = "DESeq2.pipeline1.R ../quant/featureCounts_output.txt {}".format(paramObj.ngs_sample_cond)
    subprocess.call(cmd, shell=True)

def dasOnlySR():
    pass

def processReadsFromNGS(paramObj=None, refParams=None):
    print str(datetime.datetime.now()) + " Start RNA-seq reads processing..."
    prevDir = os.getcwd()
    if os.path.exists("RNA-seq"):
        shutil.rmtree("RNA-seq")
    cmd = "mkdir -p RNA-seq/{indexDir,alignment,reassembly,quant,deg,das}"
    subprocess.call(cmd, shell=True)
    os.chdir("RNA-seq")
    gtfPrefix = os.path.splitext(refParams.ref_gff)[0]
    if refParams.ref_gtf:
        refAnnoGTF = refParams.ref_gtf
    elif refParams.ref_gff:
        refAnnoGTF = gtfPrefix + ".gtf"
        cmd = "gffread -T -o {} {}".format(refAnnoGTF, refParams.ref_genome)
        subprocess.call(cmd, shell=True)
    else:
        raise Exception("No annotaion file provided, please input the right annotion file, e.g. xx.gtf or xx.gff")
    # if defaultCfg.refParams.ref_gtf:
    #     refAnnoGTF = defaultCfg.refParams.ref_gtf
    # elif defaultCfg.refParams.ref_gff:
    #     refAnnoGTF = os.path.splitext(defaultCfg.refParams.ref_gff)[0] + ".gtf"
    #     cmd = "gffread -T -o $refAnnoGTF $refAnno"

    # transFasta = refAnnoGTF + ".fa"
    # cmd = "gffread {} -g {} -w {}".format(refParams.ref_gtf, refParams.ref_genome, transFasta)
    # subprocess.call(cmd, shell=True)
    # cmd = "salmon index -t {} -i salmon_index -type quasi -k 31".format(transFasta)
    # subprocess.call(cmd, shell=True)

    cmd = "hisat2_extract_splice_sites.py {} >{}.ss".format(refAnnoGTF, gtfPrefix)
    subprocess.call(cmd, shell=True)
    cmd = "hisat2_extract_exons.py {} >{}.exon".format(refAnnoGTF, gtfPrefix)
    subprocess.call(cmd, shell=True)
    cmd = "hisat2_build -p {} --ss {}.ss --exon {}.exon {} indexDir/{}".format(thread, gtfPrefix, gtfPrefix,
                                                                               refParams.ref_genome, gtfPrefix)
    subprocess.call(cmd, shell=True)
    srFileNameDict = getSRfileNames(paramObj.ngs_dir)
    if len(srFileNameDict) > 1:
        os.chdir("alignment")
        mergeList = open("../reassembly/mergedList.txt", "w")
        bamList = []
        for srName in srFileNameDict:
            leftRead = srFileNameDict[srName]["r1"]
            rightRead = srFileNameDict[srName]["r2"]
            hisat2Para = "-x indexDir/{} -1 {} -2 {} --dta -p {} -S {}.sam".format(gtfPrefix, leftRead, rightRead, thread,
                                                                                   srName)
            cmd = "hisat2 {}".format(hisat2Para)
            subprocess.call(cmd, shell=True)
            cmd = "samtools sort -@ {} -o {}.sorted.bam {}.sam".format(thread, srName, srName)
            subprocess.call(cmd, shell=True)
            cmd = "stringtie {}.sorted.bam -G {} -o {}.gtf -p {}".format(srName, refAnnoGTF, srName, thread)
            subprocess.call(cmd, shell=True)
            # cmd = '''echo {}.gtf >> mergedList.txt'''.format(srName)
            # subprocess.call(cmd, shell=True)
            cmd = "samtools index {}.sorted.bam".format(srName)
            subprocess.call(cmd, shell=True)
            bamList.append("{}/{}.sorted.bam".format(os.getcwd(), srName))
            print >>mergeList, "{}.gtf".format(srName)
            #
            # cmd = "salmon quant -i salmon_index -p {} --libType A -1 {} -2 {} -o {}".format(thread, leftRead, rightRead, srName)
            # subprocess.call(cmd, shell=True)
        mergeList.close()
        os.chdir("../quant")
        bamStr = " ".join(bamList)
        cmd = "featureCounts -a {} -o featureCounts_output.txt -G {} -J {}".format(refAnnoGTF, refParams.ref_genome,
                                                                                bamStr)
        subprocess.call(cmd, shell=True)
        os.chdir("../reassembly")
        cmd = "stringtie --merge -p {} -G {} -o stringtie_merged.gtf {}".format(thread, refAnnoGTF, mergeList)
        subprocess.call(cmd, shell=True)
        os.chdir("../deg")
        # degOnlySR(bamList, paramObj=paramObj)
        os.chdir("../das")
        # dasOnlySR()
    elif len(srFileNameDict) == 1:
        pass
    else:
        raise Exception("You should put some SGS reads in {}".format(paramObj.ngs_dir))
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " RNA-seq read processing end!"

    # quantOnlySR()
    # degOnlySR()
    # detOnlySR()
    # dasOnlySR()

def processRawReads(paramObj=None, refParams=None, ccsParams=None, mecatParams=None):
    # projectName, group, tgs_dir, ngs_dir, tgs_plat, pac_plat, ngs_pair = \
    #     paramObj.project, paramObj.group, paramObj.tgs_dir, paramObj.ngs_dir, \
    #     paramObj.tgs_plat, paramObj.pac_plat, paramObj.ngs_pair
    tgs_dir, tgs_plat, ngs_dir, projectName = paramObj.tgs_dir, paramObj.tgs_plat, paramObj.ngs_dir, paramObj.project

    # outDir = refParmas.out_dir
    # workDir = os.path.join(curDir, outDir, projectName)
    # if os.path.exists(outDir):
    #     shutil.rmtree(outDir)
    # os.makedirs(outDir)
    # os.chdir(workDir)

    paramObj.tgs_dir = os.path.join(originDir, paramObj.tgs_dir)
    if paramObj.ngs_dir:
        paramObj.ngs_dir = os.path.join(originDir, paramObj.ngs_dir)

    # if tgs_plat and ngs_dir == None:
    #     if tgs_plat.lower() == "pacbio":
    #         raise Exception("PacBio reads only can't be quantified!")
    #     elif tgs_plat.lower() == "nanopore":
    #         print "Nanopore can be quantified"
    # elif tgs_plat and ngs_dir:
    #     if tgs_plat.lower() == "pacbio":
    #         print "It's PACBIO hybrid strategy"
    #     elif tgs_plat.lower() == "nanopore":
    #         print "It's NANOPORE hybrid strategy"
    # elif tgs_plat == None and ngs_dir:
    #     print "It's NGS strategy"

    if tgs_dir != None:
        if tgs_plat.lower() == "pacbio":
            inputBam, seqNeedCorrect, subreads = processReadsFromPacbio(paramObj=paramObj, ccsParams=ccsParams)
            # inputBam, seqNeedCorrect, subreads = "CCS.bam", "paTrimmed.fa", "input.subreadset.xml"
            correctWithDifferentTools(seqNeedCorrect, paramObj=paramObj, refParams=refParams, mecatParams=mecatParams,
                                      inputBam=inputBam, subreads=subreads, ngs_dir=ngs_dir, tgs_plat="pacbio")
        elif tgs_plat.lower() == "nanopore":
            seqNeedCorrect = "nanopore.fq"
            correctWithDifferentTools(seqNeedCorrect, paramObj=paramObj, refParams=refParams, mecatParams=mecatParams,
                                      ngs_dir=ngs_dir, tgs_plat="nanopore")
        else:
            raise Exception("Please input the right TGS sequencing strategy, PacBio or Nanopore!")

    elif ngs_dir != None:
        processReadsFromNGS(paramObj=paramObj, refParams=refParams)
    else:
        raise Exception("Please input the data you want to be analyzed. "
                        "By setting the value of 'tgs_dir' or 'ngs_dir' in research sections")

def contamination(assignBed, refParams=None):
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
    # exomeLength = float(os.popen(binPath + "/gpeFeature.pl -e " + refGPE + "| awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}'").read().strip())
    # intromeLength = float(os.popen(binPath + "/gpeFeature.pl -i " + refGPE + "| awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}'").read().strip())
    genomeLength = 0
    with open(refSize) as f:
        for line in f:
            genomeLength += int(line.strip("\n").split("\t")[1])
    intergenomeLength = genomeLength - exomeLength - intromeLength
    exomeRPKM = float(readFeatureDict["exonic"])/(exomeLength/1000)/(float(readCount)/1000000)
    intromeRPKM = float(readFeatureDict["intronic"])/(intromeLength/1000)/(float(readCount)/1000000)
    intergenomeRPKM = float(readFeatureDict["intergenic"])/(intergenomeLength/1000)/(float(readCount)/1000000)
    log = open("contamination.log", "w")
    print >>log, "Exome\t{}".format(exomeRPKM)
    print >>log, "Introme\t{}".format(intromeRPKM)
    print >>log, "Intergenome\t{}".format(intergenomeRPKM)
    log.close()

def saturation(bedFile, refParams=None):
    cmd = "readsAssigner.pl -g {} {} >assign.multLine.bed12+".format(refParams.ref_gpe, bedFile)
    subprocess.call(cmd, shell=True)
    assignBedDict = {}
    with open("assign.multLine.bed12+") as f:
        for line in f:
            infoList = line.strip().split("\t")
            assignBedDict[infoList[3]] = infoList[14:17]
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
    for i in xrange(0, halfOfTotal, step1):
        random.shuffle(bedLines)
        headIreads = bedLines[0:i]
        shufOut = open("shuf.bed12+", "w")
        for z in headIreads:
            print >>shufOut, z.strip()
        shufOut.close()
        shufReads = [b.strip().split("\t")[3] for b in headIreads]
        transReads = set()
        geneReads = set()
        for j in shufReads:
            if j in assignBedDict and assignBedDict[j][0] == "E":
                transReads.add(assignBedDict[j][1])
                geneReads.add(assignBedDict[j][2])
        transDiscoverRate = len(transReads)/float(transN)
        print >>transDis, "{}\t{}".format(i, transDiscoverRate)
        geneDiscoverRate = len(geneReads)/float(geneN)
        print >>geneDis, "{}\t{}".format(i, geneDiscoverRate)
        cmd = '''(printf "{}\t"; exonDiscoverRate.pl -g {} shuf.bed12+) >>exon.discover'''.format(i, refParams.ref_gpe)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
    for i in xrange(halfOfTotal, totalReadN, step2):
        random.shuffle(bedLines)
        headIreads = bedLines[0:i]
        shufOut = open("shuf.bed12+", "w")
        for z in headIreads:
            print >>shufOut, z.strip()
        shufOut.close()
        shufReads = [b.strip().split("\t")[3] for b in headIreads]
        transReads = set()
        geneReads = set()
        for j in shufReads:
            if j in assignBedDict and assignBedDict[j][0] == "E":
                transReads.add(assignBedDict[j][1])
                geneReads.add(assignBedDict[j][2])
        transDiscoverRate = len(transReads)/float(transN)
        print >>transDis, "{}\t{}".format(i, transDiscoverRate)
        geneDiscoverRate = len(geneReads)/float(geneN)
        print >>geneDis, "{}\t{}".format(i, geneDiscoverRate)
        cmd = '''(printf "{}\t"; exonDiscoverRate.pl -g {} shuf.bed12+) >>exon.discover'''.format(i, refParams.ref_gpe)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
    # os.remove("shuf.bed12+")
    transDis.close()
    geneDis.close()
    exonDis.close()

def postMappingEvaluation(refParams=None):
    print str(datetime.datetime.now()) + " Post-mapping evaluation"
    prevDir = os.getcwd()
    if os.path.exists("postEvaluation/length"):
        shutil.rmtree("postEvaluation")
    os.makedirs("postEvaluation/length")
    os.chdir("postEvaluation")
    # Coverage across Genes with Differential Length
    os.symlink("../collapse/processed.ignore_id_removed.bed12+", "processed.ignore_id_removed.bed12+")
    os.symlink("../mapping/processed.fa", "processed.fa")
    # cmd = "ln -sf ../mapping/processed.ignore_id_removed.bed12+ && ln -sf ../collapse/sqanti_artifact_removed.ignore_id_removed.fa"
    cmd = "faCount processed.fa > processed.faCount"
    subprocess.call(cmd, shell=True)
    readLengthOut = open("length/Read.lst", "w")
    with open("processed.faCount") as f:
        for line in f:
            if line.startswith("#") or "total" in line: continue
            print >>readLengthOut, line.strip().split("\t")[1]
    readLengthOut.close()
    # cmd = '''grep -v '^#' sqanti_artifact_removed.ignore_id_removed.faCount | grep -v total | cut -f2 >length/Read.lst'''
    cmd = "gpe2bed.pl {} | bedLength.pl | cut -f13 >length/ReferenceGene.lst".format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True)
    cmd = '''cut -f1-12 processed.ignore_id_removed.bed12+ | bedLength.pl | cut -f13 >length/processed.ignore_id_removed.lst'''
    subprocess.call(cmd, shell=True)
    cmd = '''distrCurves.R -x1=0 -x2=10000 -d -x='Binned Length (limited in 0-10000)' -w=15 length/*.lst -b=150 -p=LengthDistribution.curve.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "boxes.R -ng -no length/*.lst -p=LengthDistribution.box.pdf 2>/dev/null"
    subprocess.call(cmd, shell=True)
    cmd = "geneCoverage.pl -g {} processed.ignore_id_removed.bed12+ >coverage.log".format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True)
    cmd = '''awk '$3<20000' coverage.log | select.pl --index 3,2 | binBox.R -b=10 -w=15 -m='Coverage across Genes' -x='RefSeq Gene Length (<20000)' -y=Coverage -p=coverage.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    # Distance to TSS/TTS across Reads with Different Length
    cmd = '''distanceToEnd.pl -g {} -f processed.faCount processed.ignore_id_removed.bed12+ | awk '$2!="NA"' >distance2Gene.log'''.format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''select.pl --index 4,2 distance2Gene.log | binBox.R -ng -no -b=5 -w=15 -m='Distance to Gene TSS' -x='Binned Read Length' -y=Distance -p=distance2TSS.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''select.pl --index 4,3 distance2Gene.log | binBox.R -ng -no -b=5 -w=15 -m='Distance to Gene TTS' -x='Binned Read Length' -y=Distance -p=distance2TTS.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "cut -f2 distance2Gene.log >DistanceToTSS.dis && cut -f3 distance2Gene.log >DistanceToTTS.dis"
    subprocess.call(cmd, shell=True)
    cmd = '''distrCurve.R <DistanceToTTS.dis -d -x='Distance (limited in -500 to 500)' -y='Density of Reads' -m='Distance to TTS' -x1=-500 -x2=500 -b=10 -p=distance2EndOfTTS.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    cmd = "rm *.dis"
    subprocess.call(cmd, shell=True)

    # Expression in different genome region
    cmd = "readsAssigner.pl -g {} -s processed.ignore_id_removed.bed12+ >assign.bed12+".format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True)
    contamination("assign.bed12+", refParams=refParams)
#     cmd = '''function genomeAbundance(){
# 	exonicReadN=$(awk '$15=="E"' $1 | wc -l)
# 	intronicReadN=$(awk '$15=="I"' $1 | wc -l)
# 	intergenicReadN=$(awk '$15=="IG"' $1 | wc -l)
# 	totalReadN=$(wc -l <$1)
# 	exomeLength=$(gpeFeature.pl -e $2 | awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}')
# 	intromeLength=$(gpeFeature.pl -i $2 | awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}')
# 	genomeLength=$(awk 'BEGIN{SUM=0}{SUM+=$2}END{print SUM}' $3)
# 	genomeLength=$(($genomeLength * 2))	### why multiple 2
# 	intergenomeLength=$(($genomeLength - $exomeLength - $intromeLength))
# 	exomeRPKM=$(echo "scale=3;$exonicReadN / ($exomeLength/1000) / ($totalReadN / 10^6)" | bc)
# 	intromeRPKM=$(echo "scale=3;$intronicReadN / ($intromeLength/1000) / ($totalReadN / 10^6)" | bc)
# 	intergenomeRPKM=$(echo "scale=3;$intergenicReadN / ($intergenomeLength/1000) / ($totalReadN / 10^6)" | bc)
# 	echo -e "Exome\t$exomeRPKM"
# 	echo -e "Introme\t$intromeRPKM"
# 	echo -e "Intergenome\t$intergenomeRPKM"
# }'''
#     cmd = "genomeAbundance assign.bed12+ ../../raw.gpe $chrSize >genomeAbundance.log"
    cmd = '''bar.R <contamination.log -anno -c=darkgreen -f=white -m=Contamination -x='Genomic Context' -y=RPKM -p=contamination.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)

    # Saturation of Sequencing
    saturation("processed.ignore_id_removed.bed12+", refParams=refParams)
#     cmd = "totalReadN=$(wc -l <processed.ignore_id_removed.bed12+)"
#     cmd = '''function saturation(){
# 	readsAssigner.pl -g $1 $2 >assign.multLine.bed12+
# 	transN=$(wc -l <$1)
# 	geneN=$(cut -f12 $1 | sort -u | wc -l)
# 	step1=$(($totalReadN/80))
# 	step2=$(($totalReadN/40))
# 	halfOfTotal=$(($totalReadN/2))
# 	printf '' >trans.discover
# 	printf '' >gene.discover
# 	printf '' >exon.discover
# 	for(( i=0; i<=$halfOfTotal; i=$(($i+$step1)) ));do
# 		shuf $2 | head -n $i >shuf.bed12+
# 		transReads=$(filter.pl -o shuf.bed12+ -1 4 -2 4 -m i assign.multLine.bed12+ | awk '$15=="E"{print $16}' | sort -u | wc -l)
# 		geneReads=$(filter.pl -o shuf.bed12+ -1 4 -2 4 -m i assign.multLine.bed12+ | awk '$15=="E"{print $17}' | sort -u | wc -l)
# 		discoverRate=$(echo "scale=3;$transReads/$transN" | bc)
# 		echo -e "$i\t$discoverRate" >>trans.discover
# 		discoverRate=$(echo "scale=3;$geneReads/$geneN" | bc)
# 		echo -e "$i\t$discoverRate" >>gene.discover
# 		printf "$i\t" >>exon.discover
# 		exonDiscoverRate.pl -g $1 shuf.bed12+ >>exon.discover
# 	done
# 	for((; i<=$totalReadN; i=$(($i+$step2)) ));do
# 		shuf $2 | head -n $i >shuf.bed12+
# 		transReads=$(filter.pl -o shuf.bed12+ -1 4 -2 4 -m i assign.multLine.bed12+ | awk '$15=="E"{print $16}' | sort -u | wc -l)
# 		geneReads=$(filter.pl -o shuf.bed12+ -1 4 -2 4 -m i assign.multLine.bed12+ | awk '$15=="E"{print $17}' | sort -u | wc -l)
# 		discoverRate=$(echo "scale=3;$transReads/$transN" | bc)
# 		echo -e "$i\t$discoverRate" >>trans.discover
# 		discoverRate=$(echo "scale=3;$geneReads/$geneN" | bc)
# 		echo -e "$i\t$discoverRate" >>gene.discover
# 		printf "$i\t" >>exon.discover
# 		exonDiscoverRate.pl -g $1 shuf.bed12+ >>exon.discover
# 	done
# 	rm -f shuf.bed12+
# }'''
#     cmd = "saturation ../../raw.gpe processed.ignore_id_removed.bed12+"
    cmd = '''lines.R *.discover -m='Saturation of Sequencing' -x='Reads Count' -y='Discovery Rate' -w=15 -lgPos=top -p=discoveryRate.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)

    # Sequencing Error Rate
    errorRate = open("errorRate.log", "w")
    with open("processed.ignore_id_removed.bed12+") as f:
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            print >>errorRate, "{}\t{}".format(lineInfo[3], 100-float(lineInfo[12])*float(lineInfo[13])/100)
    errorRate.close()
    cmd = '''cut -f2 errorRate.log | distrCurve.R - -d -m='Rate of Mismatch & Indel for Mapping' -x=Rate -y='Fraction of Reads' -p=errorRate.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    # cmd = '''awk 'BEGIN{OFS="\t"}{print $4,100-$13*$14/100}' processed.ignore_id_removed.bed12+ | tee errorRate.log | cut -f2 | distrCurve.R -d -m='Rate of Mismatch & Indel for Mapping' -x=Rate -y='Fraction of Reads' -p=errorRate.pdf 2>/dev/null'''
    # Read Count per Gene
    cmd = "geneRPKM.pl -g {} processed.ignore_id_removed.bed12+ | tee RPKM.bed6+ | cut -f7 >readCount.lst".format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True)
    cmd = '''distrCurve.R <readCount.lst -m='Gene Count at Read Count' -x='Read Count' -y='Gene Count' -xl=10 -yl=10 -ng -p=readCountPerGene.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)

    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Post-mapping evaluation done!"

def paAnalysis(refParams=None):
    print str(datetime.datetime.now()) + " PA Analysis..."
    prevDir = os.getcwd()
    if os.path.exists("PA"):
        shutil.rmtree("PA")
    os.makedirs("PA")
    os.chdir("PA")
    os.symlink("../collapse/processed.ignore_id_removed.bed12+", "processed.ignore_id_removed.bed12+")
    os.symlink("../mapping/processed.bed12+", "processed.bed12+")
    with open("processed.ignore_id_removed.bed12+") as f:
        cleaveList = []
        count = 0
        cleaveOut = open("cleavage.bed6", "w")
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            chrom, start, end, strand = lineInfo[0], int(lineInfo[1]), int(lineInfo[2]), lineInfo[5]
            if strand == "+":
                keyStr = "{}:{}-{}:{}".format(chrom, end-1, end, strand)
                if keyStr not in cleaveList:
                    count += 1
                    cleaveList.append(keyStr)
                    print >>cleaveOut, "\t".join(map(str, [chrom, end - 1, end, "ClevageSite{}".format(count), 0, strand]))
            elif strand == "-":
                keyStr = "{}:{}-{}:{}".format(chrom, start, start + 1, strand)
                if keyStr not in cleaveList:
                    count += 1
                    cleaveList.append(keyStr)
                    print >>cleaveOut, "\t".join(map(str, [chrom, start, start + 1, "ClevageSite{}".format(count), 0, strand]))
        cleaveOut.close()

    # cmd = '''awk 'BEGIN{FS=OFS="\t"}$6=="+"{$2=$3-1}$6=="-"{$3=$2+1}{print $1,$2,$3,$6}' processed.ignore_id_removed.bed12+ | sort -u | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"ClevageSite"NR,0,$4}' >cleavage.bed6'''
    cmd = '''bedClosest.pl cleavage.bed6 | tee adjacentCleavage.tsv | cut -f13 | distrCurve.R -xl=2 -d -x='Binned Distance between Adjacent Cleavage Site' -p=adjacentCleavageDistr.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    try:
        os.remove("paCluster.count")
    except OSError:
        pass
    disList = range(1, 100) + range(100, 10000, 10) + range(10000, 100000, 1000)
    paClusterCount = open("paCluster.count", "w")
    for i in disList:
        count = int(os.popen(binPath+"/paCluster.pl -d {} -m mode -w 3 processed.bed12+ | wc -l".format(i)).read().strip())
        print >>paClusterCount, "{}\t{}".format(i, count)
    paClusterCount.close()
#     cmd = '''for i in $(seq 5 1 100) $(seq 100 10 10000) $(seq 10000 1000 100000);do
#     echo -e $i"\t"$(paCluster.pl -d $i -m mode -w 3 processed.ignore_id_removed.bed12+ | wc -l) >>paCluster.count
# done'''
#     cmd = "line.R <paCluster.count -x=Distance -y=Count -p=paClusterCount.pdf 2>/dev/null"
#     subprocess.call(cmd, shell=True)
    paDistance = 30
    cmd = '''paCluster.pl -d 30 -m mode -w 3 processed.ignore_id_removed.bed12+ | tee paCluster.bed8+ | awk '$5>1{print $3-$2}' | distrCurve.R -d -x='Cluster Size (limited in 1-100)' -y='Density' -m='Distribution of Cluster Size' -x1=0 -x2=100 -b=1 -p=paClusterSize.pdf'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    paBed6 = open("PA.bed6+", "w")
    with open("paCluster.bed8+") as f:
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            print >>paBed6, "\t".join(map(str, [lineInfo[0], lineInfo[6], lineInfo[7]] + lineInfo[3:6] + lineInfo[8:]))
    paBed6.close()
    # cmd = '''awk 'BEGIN{FS=OFS="\t"}{$2=$7;$3=$8;print}' paCluster.bed8+ | cut -f1-6,9- >PA.bed6+'''
    cmd = "PAClassbyRead.pl -a ../postEvaluation/assign.bed12+ <(cut -f1-8 paCluster.bed8+) >paCluster.type.bed8+ 2>singleExonReadWithExonInMEread.bed12+"
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    deSingleExon = open("paCluster.deSingleExonRead.bed8+", "w")
    singleExon = open("paCluster.singleExonRead.bed8+", "w")
    with open("paCluster.type.bed8+") as f:
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            if lineInfo[8] != "SE":
                print >>deSingleExon, line.strip("\n")
            elif lineInfo[8] == "SE":
                print >>singleExon, line.strip("\n")
    deSingleExon.close()
    singleExon.close()
    # cmd = '''awk '$9!="SE"' paCluster.type.bed8+ >paCluster.deSingleExonRead.bed8+ && awk '$9=="SE"' paCluster.type.bed8+ >paCluster.singleExonRead.bed8+'''
    try:
        shutil.rmtree("motif")
    except OSError:
        pass
    os.makedirs("motif")
    up = 100
    down = 100
    for seq in ["A", "T", "C", "G"]:
        cmd = "awk '$5>1' PA.bed6+ | motifAroundPA.pl -u {} -d {} -m {} -b {} >motif/{}.nucleotide 2>/dev/null".format(up, down, seq, refParams.ref2bit, seq)
        subprocess.call(cmd, shell=True)
#     cmd = '''for seq in A T C G;do
#     awk '$5>1' PA.bed6+ | motifAroundPA.pl -u $up -d $down -m $seq -b $twoBit >motif/$seq.nucleotide 2>/dev/null
# done'''
    cmd = '''lines.R -w=15 -y=Frequency -x='Distance Relative to PA' -m='Distribution of Nucleotide Frequency around PA' -p=motif/nucleotide.pdf motif/*.nucleotide 2>/dev/null'''
    subprocess.call(cmd, shell=True)

    up2 = 50
    down2 = 0
    for seq in ["AATAAA", "AAATAA", "ATAAAA", "ATTAAA", "ATAAAT", "TAATAA", "ATAAAG", "AAAATA", "CAATAA", "ATAAAC", "AAAAAA", "AAAAAG"]:
        cmd = "awk '$5>1' PA.bed6+ | motifAroundPA.pl -u {} -d {} -m {} -b {} >motif/{}.PAS 2>/dev/null".format(up2, down2, seq, refParams.ref2bit, seq)
        subprocess.call(cmd, shell=True)
#     cmd = '''for seq in AATAAA AAATAA ATAAAA ATTAAA ATAAAT TAATAA ATAAAG AAAATA CAATAA ATAAAC AAAAAA AAAAAG;do
#     awk '$5>1' PA.bed6+ | motifAroundPA.pl -u $up2 -d $down2 -m $seq -b $twoBit >motif/$seq.PAS 2>/dev/null
# done'''
    cmd = '''lines.R -p=motif/PAS.1.pdf motif/{AATAAA,AAATAA,ATAAAA,ATTAAA,ATAAAT,TAATAA}.PAS -x1=-50 -x2=0 -w=15 2>/dev/null && lines.R -p=motif/PAS.2.pdf motif/{ATAAAG,AAAATA,CAATAA,ATAAAC,AAAAAA,AAAAAG}.PAS -x1=-50 -x2=0 -w=15 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''paste motif/ATAAAC.PAS <(cut -f2 motif/ATAAAG.PAS) <(cut -f2 motif/CAATAA.PAS) <(cut -f2 motif/AAAAAA.PAS) <(cut -f2 motif/AAAAAG.PAS) <(cut -f2 motif/AAAATA.PAS) | awk '{print $1"\t"$2+$3+$4+$5+$6+$7}' >motif/Other.PAS'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''lines.R -p=motif/PAS.pdf motif/{Other,AATAAA,AAATAA,ATAAAA,ATTAAA,ATAAAT,TAATAA}.PAS -x1=-50 -x2=0 -w=15 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " PA Analysis done!"

def getColFromFile(myFile, col, sep=None):
    with open(myFile) as f:
        myDict = {}
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            if sep:
                myDict.fromkeys(lineInfo[col - 1].split(sep), "")
            else:
                myDict[lineInfo[col - 1]] = ""
        return myDict

def getAnnoASList(inFile, outFile, PA=False):
    with open(inFile) as f:
        out = open(outFile, "w")
        if not PA:
            for line in f:
                if "Novel" in line: continue
                lineInfo = line.strip("\n").split("\t")
                print >>out, "{}:{}:{}".format(lineInfo[0], lineInfo[3], lineInfo[5])
        else:
            for line in f:
                lineInfo = line.strip("\n").split("\t")
                paPos = lineInfo[3].split(",")
                if len(paPos) > 1:
                    print >>out, ":".join(lineInfo[0:4])
        out.close()

def filter(originFile=None, targetFile=None, originField=1, targetField=1, mode="i", outFile=None):
    originList, targetList = [], []
    if isinstance(originFile, list):
        originList = originFile
    elif os.path.isfile(originFile):
        with open(originFile) as f:
            originList = [line.strip("\n").split("\t")[originField - 1] for line in f.readlines()]

    returnList = []
    outFile = outFile if outFile != None else "filtered.txt"
    out = open(outFile, "w")
    if isinstance(targetFile, list):
        targetList = targetFile
        if mode == "i":
            for i in targetList:
                if i in originList:
                    print >>out, i
                    returnList.append(i)
        else:
            for i in targetList:
                if i not in originList:
                    print >>out, i
                    returnList.append(i)
    elif isinstance(targetFile, dict):
        if mode == "i":
            for i in targetFile:
                if i in originList:
                    print >>out, targetFile[i]
                    returnList.append(targetFile[i])
        else:
            for i in targetFile:
                if i not in originList:
                    print >>out, targetFile[i]
                    returnList.append(targetFile[i])
    elif os.path.isfile(targetFile):
        with open(targetFile) as f:
            if mode == "i":
                for line in f:
                    lineInfo = line.strip("\n").split("\t")
                    if lineInfo[targetField-1] in originList:
                        print >>out, line.strip("\n")
                        returnList.append(line.strip("\n"))
            else:
                for line in f:
                    lineInfo = line.strip("\n").split("\t")
                    if lineInfo[targetField-1] not in originList:
                        print >>out, line.strip("\n")
                        returnList.append(line.strip("\n"))
    out.close()
    return returnList

def getASstatistics(asType="IR", asFile=None, annoFile=None, novelFile=None, outFile=None):
    out = open(outFile, "w")
    if asType == "PA":
        print >> out, "#Chr\tStrand\tKnown or Novel\tGene\tPA Sites"
        with open(annoFile) as f1:
            for line in f1:
                lineInfo = line.strip("\n").split("\t")
                if len(lineInfo[3].split(",")) > 1:
                    print >>out, "\t".join([lineInfo[0], lineInfo[1], "Known", lineInfo[2], lineInfo[3]])
        with open(novelFile) as f2:
            for line in f2:
                lineInfo = line.strip("\n").split("\t")
                if len(lineInfo[3].split(",")) > 1:
                    print >>out, "\t".join([lineInfo[0], lineInfo[1], "Novel", lineInfo[2], lineInfo[3]])
    else:
        asDict = {}
        with open(asFile) as f:
            for line in f:
                lineInfo = line.strip("\n").split("\t")
                key = ":".join([lineInfo[0], lineInfo[3], lineInfo[5]])
                asDict[key] = "\t".join([key, line.strip("\n")])
        annoList = filter(originFile=annoFile, targetFile=asDict)
        novelList = filter(originFile=novelFile, targetFile=asDict)
        os.remove("filtered.txt")
        geneCol = 7
        if asType == "IR":
            print >> out, "##AS ID is composed of Gene:Retained intron start-Retained intron end"
        elif asType == "SE":
            geneCol = 13
            print >> out, "##AS ID is composed of Gene:Left flanking constitutive exon end@Alternative exons locus@Right flanking constitutive exon start"
            print >> out, "##Alternative exons locus is composed of Alternative exon1 start-Alternative exon1 end[;Alternative exon2 start-Alternative exon2 end[;Alternative exon3 start-Alternative exon3 end...]"
        elif asType == "A3SS":
            print >> out, "##AS ID is composed of Gene:Alternative 5' splicing region start-Alternative 5' splicing region end"
        elif asType == "A5SS":
            print >> out, "##AS ID is composed of Gene:Alternative 3' splicing region start-Alternative 3' splicing region end"
        print >> out, "#Chr\tStrand\tKnown or Novel\tAS ID\tGene"

        for i in annoList:
            tmpList = i.strip().split("\t")
            print >> out, "\t".join([tmpList[1], tmpList[6], "Known", tmpList[4], tmpList[geneCol]])
        for i in novelList:
            tmpList = i.strip().split("\t")
            print >> out, "\t".join([tmpList[1], tmpList[6], "Novel", tmpList[4], tmpList[geneCol]])

    out.close()

def identifyASE(refParams=None, paramObj=None, project2as=None):
    # Identify alternative splicing events
    print str(datetime.datetime.now()) + " Alternative splicing events identifying..."
    prevDir = os.getcwd()
    if os.path.exists("ASE"):
        shutil.rmtree("ASE")
    os.makedirs("ASE")
    os.chdir("ASE")
    cmd = "readGroup.pl -g {} ../postEvaluation/assign.bed12+ >unambi.bed12+ 2>ambiguous.bed12+".format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True)
    with open("unambi.bed12+") as f:
        deSingleExonRead = open("deSingleExonRead.bed12+", "w")
        readGroup = open("readGrouped.bed12+", "w")
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            meDict = getColFromFile("../PA/singleExonReadWithExonInMEread.bed12+", 4, sep=",")
            seDict = getColFromFile("../PA/paCluster.singleExonRead.bed8+", 4, sep=",")
            if lineInfo[3] not in meDict and lineInfo[3] not in seDict:
                print >>deSingleExonRead, "\t".join(lineInfo[0:12]+lineInfo[20:21])
            print >>readGroup, "\t".join(lineInfo[0:12]+lineInfo[20:21])
        deSingleExonRead.close()
        readGroup.close()
    # cmd = '''cut -f1-12,21 unambi.bed12+ | tee readGrouped.bed12+ | filter.pl -o <(cut -f4 ../PA/singleExonReadWithExonInMEread.bed12+; cut -f4 ../PA/paCluster.singleExonRead.bed8+ | tr ',' "\n") -2 4 >deSingleExonRead.bed12+'''
    cmd = "3endRevise.pl -p ../PA/PA.bed6+ deSingleExonRead.bed12+ | tee deSingleExonRead.3endRevised.bed12+ | paGroup.pl >paGrouped.tsv 2>paGrouped.bed6"
    subprocess.call(cmd, shell=True)
    os.makedirs("PB")
    findASmain(asType="IR", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/IR.bed6+")
    # cmd = "getSingleIr.pl -g {} deSingleExonRead.bed12+ >PB/IR.bed6+".format(refParams.ref_gpe)
    # subprocess.call(cmd, shell=True)
    geneWithIR = len(getColFromFile("PB/IR.bed6+", 7))
    deSingleExonReadCount = len(getColFromFile("deSingleExonRead.bed12+", 13))
    geneWithIrRatio = round(float(geneWithIR)/deSingleExonReadCount, 4)
    # cmd = "geneWithIR=$(cut -f7 PB/IR.bed6+ | sort -u | wc -l)"
    # cmd = "deSingleExonReadCount=$(cut -f13 deSingleExonRead.bed12+ | sort -u | wc -l)"
    # cmd = '''geneWithIrRatio=$(awk 'BEGIN{printf "%.4f", '$geneWithIR/$deSingleExonReadCount'}')'''
    # path2findAS = "/home/xufeng/xufeng/iso-seq/test/workflow_scripts_test/test_20190808/scripts"
    # cmd = "python {}/findASmain.py -g {} -a {} -s {} -o {}".format(path2findAS, refParams.ref_gpe, "SE", "deSingleExonRead.bed12+",
    #                                                      "PB/SE.bed12+")
    # subprocess.call(cmd, shell=True)
    # cmd = "python {}/findASmain.py -g {} -a {} -s {} -o {}".format(path2findAS, refParams.ref_gpe, "A3SS", "deSingleExonRead.bed12+",
    #                                                      "PB/A3SS.bed6+")
    # subprocess.call(cmd, shell=True)
    # cmd = "python {}/findASmain.py -g {} -a {} -s {} -o {}".format(path2findAS, refParams.ref_gpe, "A5SS", "deSingleExonRead.bed12+",
    #                                                      "PB/A5SS.bed6+")
    # subprocess.call(cmd, shell=True)
    findASmain(asType="SE", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/SE.bed12+")
    findASmain(asType="A5SS", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/A5SS.bed6+")
    findASmain(asType="A3SS", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/A3SS.bed6+")
    project2as[paramObj.project] = os.path.join(os.getcwd(), "deSingleExonRead.bed12+")
    # cmd = "getSE.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/SE.bed12+ && getA5SS.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/A5SS.bed6+ && getA3SS.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/A3SS.bed6+"

    cmd = '''awk '$9>1 && $11>1 && $5>=100 && $5<=900{print $3-$2}' PB/A5SS.bed6+ | hist.R -p=PB/A5SS.size.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    cmd = '''awk '$3-$2>=20 && $9>1 && $11>1{print $5/1000}' PB/A5SS.bed6+ | hist.R -x='Incusion Ratio' -p=PB/A5SS.InclusionRatio.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    os.symlink("SE.bed12+", "PB/SE.confident.bed12+")
    os.symlink("A5SS.bed6+", "PB/A5SS.confident.bed6+")
    os.symlink("A3SS.bed6+", "PB/A3SS.confident.bed6+")
    # cmd = "ln -sf SE.bed12+ PB/SE.confident.bed12+ && ln -sf A5SS.bed6+ PB/A5SS.confident.bed6+ && ln -sf A3SS.bed6+ PB/A3SS.confident.bed6+"
    cmd = '''
            awk '{print $3-$2}' PB/A5SS.confident.bed6+ | hist.R -p=PB/A5SS.size.pdf 2>/dev/null
            awk '$9>1 && $11>1{print $5/1000}' PB/A5SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=PB/A5SS.InclusionRatio.pdf 2>/dev/null
            awk '{print $3-$2}' PB/A3SS.confident.bed6+ | hist.R -p=PB/A3SS.size.pdf 2>/dev/null
            awk '$9>1 && $11>1{print $5/1000}' PB/A3SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=PB/A3SS.InclusionRatio.pdf 2>/dev/null
        '''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    print str(datetime.datetime.now()) + " Alternative splicing events identifying done!"

    # Alternative splicing events characterization
    print str(datetime.datetime.now()) + " ASE Characterization..."
    os.makedirs("characterization")
    os.chdir("characterization")
    getAnnoASList("../PB/SE.bed12+", "PB.SE.lst")
    getAnnoASList("../PB/A5SS.bed6+", "PB.A5SS.lst")
    getAnnoASList("../PB/A3SS.bed6+", "PB.A3SS.lst")
    getAnnoASList("../PB/IR.bed6+", "PB.IR.lst")
    getAnnoASList("../paGrouped.tsv", "PB.APA.lst", PA=True)
    # cmd = '''
    #         awk '{print $1":"$4":"$6}' ../PB/SE.bed12+              | grep -v Novel >PB.SE.lst
    #         awk '{print $1":"$4":"$6}' ../PB/A5SS.bed6+             | grep -v Novel >PB.A5SS.lst
    #         awk '{print $1":"$4":"$6}' ../PB/A3SS.bed6+             | grep -v Novel >PB.A3SS.lst
    #         awk '{print $1":"$4":"$6}' ../PB/IR.bed6+               | grep -v Novel >PB.IR.lst
    #         awk 'BEGIN{FS="\t";OFS=":"}{split($4,array,",")}length(array)>1&&$3!~/^chr/{print $1,$2,$3,$4}' ../paGrouped.tsv >PB.APA.lst
    #     '''
    cmd = '''
            venn.R PB.IR.lst -p=IR.pdf 2>/dev/null
            venn.R PB.APA.lst -p=APA.pdf 2>/dev/null
            ln -sf PB.SE.lst confident.SE.lst
            ln -sf PB.A5SS.lst confident.A5SS.lst
            ln -sf PB.A3SS.lst confident.A3SS.lst
        '''
    subprocess.call(cmd, shell=True)
    # Classify APE into annotated or novel according to reference gene model
    cmd = '''awk '{{print $0"\t"$12}}' {} | gpe2bed.pl -p >reference.bed12+'''.format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    # cmd = "python {}/findASmain.py -g {} -a {} -s {} -o {}".format(path2findAS, refParams.ref_gpe, "IR", "reference.bed12+",
    #                                                      "IR.reference.bed6+")
    # subprocess.call(cmd, shell=True)
    # cmd = "python {}/findASmain.py -g {} -a {} -s {} -o {}".format(path2findAS, refParams.ref_gpe, "SE", "reference.bed12+",
    #                                                      "SE.reference.bed6+")
    # subprocess.call(cmd, shell=True)
    # cmd = "python {}/findASmain.py -g {} -a {} -s {} -o {}".format(path2findAS, refParams.ref_gpe, "A3SS", "reference.bed12+",
    #                                                      "A3SS.reference.bed6+")
    # subprocess.call(cmd, shell=True)
    # cmd = "python {}/findASmain.py -g {} -a {} -s {} -o {}".format(path2findAS, refParams.ref_gpe, "A5SS", "reference.bed12+",
    #                                                      "A5SS.reference.bed6+")
    # subprocess.call(cmd, shell=True)
    # subprocess.call(cmd, shell=True, stderr=sys.stderr, stdout=sys.stdout)

    findASmain(asType="IR", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="IR.reference.bed6+")
    findASmain(asType="SE", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="SE.reference.bed6+")
    findASmain(asType="A3SS", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="A3SS.reference.bed6+")
    findASmain(asType="A5SS", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="A5SS.reference.bed6+")
    getAnnoASList("IR.reference.bed6+", "IR.reference.lst")
    getAnnoASList("SE.reference.bed6+", "SE.reference.lst")
    getAnnoASList("A3SS.reference.bed6+", "A3SS.reference.lst")
    getAnnoASList("A5SS.reference.bed6+", "A5SS.reference.lst")
    cmd = "paGroup.pl reference.bed12+ >paGroup.reference.tsv 2>paGroup.reference.bed6"
    subprocess.call(cmd, shell=True)
    cmd = "paCmp.pl -r paGroup.reference.tsv -a annoPA.tsv -n novelPA.tsv ../paGrouped.tsv >paGroup.novel.tsv 2>paGroup.anno.tsv"
    subprocess.call(cmd, shell=True)
    # cmd = '''
    #         getSingleIr.pl -g ../../../raw.gpe reference.bed12+ | tee IR.reference.bed6+ | awk '{print $1":"$4":"$6}' >IR.reference.lst
    #         paGroup.pl reference.bed12+ >paGroup.reference.tsv 2>paGroup.reference.bed6
    #         paCmp.pl -r paGroup.reference.tsv -a annoPA.tsv -n novelPA.tsv ../paGrouped.tsv >paGroup.novel.tsv 2>paGroup.anno.tsv
    #         getSE.pl -g ../../../raw.gpe reference.bed12+ | tee SE.reference.bed12+ | awk '{print $1":"$4":"$6}' >SE.reference.lst
    #         getA5SS.pl -g ../../../raw.gpe reference.bed12+ | tee A5SS.reference.bed6+ | awk '{print $1":"$4":"$6}' >A5SS.reference.lst
    #         getA3SS.pl -g ../../../raw.gpe reference.bed12+ | tee A3SS.reference.bed6+ | awk '{print $1":"$4":"$6}' >A3SS.reference.lst
    #     '''
    filter(originFile="IR.reference.lst", targetFile="PB.IR.lst", outFile="IR.anno.lst", mode="i")
    filter(originFile="IR.reference.lst", targetFile="PB.IR.lst", outFile="IR.novel.lst", mode="e")
    filter(originFile="SE.reference.lst", targetFile="confident.SE.lst", outFile="SE.anno.lst", mode="i")
    filter(originFile="SE.reference.lst", targetFile="confident.SE.lst", outFile="SE.novel.lst", mode="e")
    filter(originFile="A3SS.reference.lst", targetFile="confident.A3SS.lst", outFile="A3SS.anno.lst", mode="i")
    filter(originFile="A3SS.reference.lst", targetFile="confident.A3SS.lst", outFile="A3SS.novel.lst", mode="e")
    filter(originFile="A5SS.reference.lst", targetFile="confident.A5SS.lst", outFile="A5SS.anno.lst", mode="i")
    filter(originFile="A5SS.reference.lst", targetFile="confident.A5SS.lst", outFile="A5SS.novel.lst", mode="e")
    # cmd = '''
    #         filter.pl -o IR.reference.lst PB.IR.lst >IR.novel.lst
    #         filter.pl -o IR.reference.lst PB.IR.lst -m i >IR.anno.lst
    #         filter.pl -o SE.reference.lst   confident.SE.lst        >SE.novel.lst
    #         filter.pl -o SE.reference.lst   confident.SE.lst   -m i >SE.anno.lst
    #         filter.pl -o A5SS.reference.lst confident.A5SS.lst      >A5SS.novel.lst
    #         filter.pl -o A5SS.reference.lst confident.A5SS.lst -m i >A5SS.anno.lst
    #         filter.pl -o A3SS.reference.lst confident.A3SS.lst      >A3SS.novel.lst
    #         filter.pl -o A3SS.reference.lst confident.A3SS.lst -m i >A3SS.anno.lst
    #     '''
    # Get Statistic of APE
    getASstatistics(asType="IR", asFile="../PB/IR.bed6+", annoFile="IR.anno.lst", novelFile="IR.novel.lst",
                    outFile="statistic.IR.tsv")
    getASstatistics(asType="PA", annoFile="paGroup.anno.tsv", novelFile="paGroup.novel.tsv",
                    outFile="statistic.APA.tsv")
    getASstatistics(asType="SE", asFile="../PB/SE.confident.bed12+", annoFile="SE.anno.lst", novelFile="SE.novel.lst",
                    outFile="statistic.SE.tsv")
    getASstatistics(asType="A3SS", asFile="../PB/A3SS.confident.bed6+", annoFile="A3SS.anno.lst",
                    novelFile="A3SS.novel.lst", outFile="statistic.A3SS.tsv")
    getASstatistics(asType="A5SS", asFile="../PB/A5SS.confident.bed6+", annoFile="A5SS.anno.lst",
                    novelFile="A5SS.novel.lst", outFile="statistic.A5SS.tsv")
    # cmd = '''
    #         (echo -e "##AS ID is composed of Gene:Retained intron start-Retained intron end"
    #         echo -e "#Chr\tStrand\tKnown or Novel\tAS ID\tGene";
    #         awk '{print $1":"$4":"$6"\t"$0}' ../PB/IR.bed6+ | filter.pl -o IR.anno.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Known",$5,$8}'
    #         awk '{print $1":"$4":"$6"\t"$0}' ../PB/IR.bed6+ | filter.pl -o IR.novel.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Novel",$5,$8}'
    #         ) >statistic.IR.tsv
    #         (echo -e "#Chr\tStrand\tKnown or Novel\tGene\tPA Sites"
    #         awk 'BEGIN{OFS="\t"}{split($4,array,",")}length(array)>1{print $1,$2,"Known",$3,$4}' paGroup.anno.tsv
    #         awk 'BEGIN{OFS="\t"}{split($4,array,",")}length(array)>1{print $1,$2,"Novel",$3,$4}' paGroup.novel.tsv
    #         ) >statistic.APA.tsv
    #         (echo -e "##AS ID is composed of Gene:Left flanking constitutive exon end@Alternative exons locus@Right flanking constitutive exon start"
    #         echo -e "##Alternative exons locus is composed of Alternative exon1 start-Alternative exon1 end[;Alternative exon2 start-Alternative exon2 end[;Alternative exon3 start-Alternative exon3 end...]"
    #         echo -e "#Chr\tStrand\tKnown or Novel\tAS ID\tGene"
    #         awk '{print $1":"$4":"$6"\t"$0}' ../PB/SE.confident.bed12+ | filter.pl -o SE.anno.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Known",$5,$14}'
    #         awk '{print $1":"$4":"$6"\t"$0}' ../PB/SE.confident.bed12+ | filter.pl -o SE.novel.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Novel",$5,$14}'
    #         ) >statistic.SE.tsv
    #         (echo -e "##AS ID is composed of Gene:Alternative 5' splicing region start-Alternative 5' splicing region end"
    #         echo -e "#Chr\tStrand\tKnown or Novel\tAS ID\tGene";
    #         awk '{print $1":"$4":"$6"\t"$0}' ../PB/A5SS.confident.bed6+ | filter.pl -o A5SS.anno.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Known",$5,$8}'
    #         awk '{print $1":"$4":"$6"\t"$0}' ../PB/A5SS.confident.bed6+ | filter.pl -o A5SS.novel.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Novel",$5,$8}'
    #         ) >statistic.A5SS.tsv
    #         (echo -e "##AS ID is composed of Gene:Alternative 3' splicing region start-Alternative 3' splicing region end"
    #         echo -e "#Chr\tStrand\tKnown or Novel\tAS ID\tGene";
    #         awk '{print $1":"$4":"$6"\t"$0}' ../PB/A3SS.confident.bed6+ | filter.pl -o A3SS.anno.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Known",$5,$8}'
    #         awk '{print $1":"$4":"$6"\t"$0}' ../PB/A3SS.confident.bed6+ | filter.pl -o A3SS.novel.lst -m i | awk 'BEGIN{OFS="\t"}{print $2,$7,"Novel",$5,$8}'
    #         ) >statistic.A3SS.tsv
    #     '''
    # Event Decompose and Motif Evaluation
    getSpliceSite(asType="IR", asFile="PB.IR.lst", outFile="IR.splicesite")
    getSpliceSite(asType="IR", asFile="IR.anno.lst", outFile="IR.anno.splicesite")
    getSpliceSite(asType="IR", asFile="IR.novel.lst", outFile="IR.novel.splicesite")
    getSpliceSite(asType="SE", asFile="confident.SE.lst")
    getSpliceSite(asType="A3SS", asFile="confident.A3SS.lst")
    getSpliceSite(asType="A5SS", asFile="confident.A5SS.lst")
    splicesite2figure(refParams=refParams)

    # cmd = '''
    #         awk 'BEGIN{FS=":";OFS="\t"}{split($3,pos,"-");print $1,pos[1],pos[2],$4}' PB.IR.lst    >IR.splicesite
    #         awk 'BEGIN{FS=":";OFS="\t"}{split($3,pos,"-");print $1,pos[1],pos[2],$4}' IR.anno.lst  >IR.anno.splicesite
    #         awk 'BEGIN{FS=":";OFS="\t"}{split($3,pos,"-");print $1,pos[1],pos[2],$4}' IR.novel.lst >IR.novel.splicesite
    #         seDecompose.pl        confident.SE.lst   >SE.inc.splicesite   2>SE.exc.splicesite
    #         anssDecompose.pl -n 5 confident.A5SS.lst >A5SS.inc.splicesite 2>A5SS.exc.splicesite
    #         anssDecompose.pl -n 3 confident.A3SS.lst >A3SS.inc.splicesite 2>A3SS.exc.splicesite
    #         splicesite2figure
    #     '''
    # Distance of PA to TTS
    getDist2TTS(refParams=refParams)
    # cmd = '''
    #         awk '$4!~/^chr/' ../paGrouped.bed6 | awk 'BEGIN{OFS="\t"}$6=="+"{print $1,$3-1,$3,$4,$5,$6}$6=="-"{print $1,$2,$2+1,$4,$5,$6}' >pbPA.bed6
    #         bedtools closest -a <(sort -k1,1 -k2,2n pbPA.bed6) -b <(gpeFeature.pl --tts ../../../raw.gpe|sort -k1,1 -k2,2n) -s -D a | select.pl -i 13,4 | sort -u | tee pbPA2TTS.tsv | cut -f1 | box.R -ng -nJ -no -y='Distance to TTS' -p=pbPA2TTS.pdf
    #     '''
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " ASE Characterization done!"
    return project2as

def findASmain(asType="IR", gpeFile=None, sourceFile=None, outFile=None):
    gpeObj = GenePredObj(gpeFile, False)
    gpeDict = gpeObj.genePredDict
    sampleObj = Sample()
    sample = "temp"
    with open(sourceFile) as f:
        for line in f:
            readStruc = ReadLineStruc(line)
            sampleObj.update(sample, sourceFile, readStruc)
    geneDict, novelDict = {}, {}
    for g in sampleObj.sample2Gene[sample]:
        getAS(gpeDict, sampleObj.sample2Gene[sample][g], geneDict, novelDict, offset=0)
    out = open(outFile, "w")
    findAS(geneDict, outASType=asType, anno=True, out=out)
    findAS(novelDict, outASType=asType, anno=False, out=out)
    out.close()

def getDist2TTS(refParams=None):
    with open("../paGrouped.bed6") as f:
        out = open("pbPA.bed6", "w")
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            if lineInfo[5] == "+":
                print >>out, "\t".join(map(str, [lineInfo[0], int(lineInfo[2])-1, lineInfo[2]] + lineInfo[3:]))
            else:
                print >>out, "\t".join(map(str, [lineInfo[0], lineInfo[1], int(lineInfo[1])+1] + lineInfo[3:]))
        out.close()
    cmd = "bedtools closest -a <(sort -k1,1 -k2,2n pbPA.bed6) -b <(gpeFeature.pl --tts {}|sort -k1,1 -k2,2n) -s -D a | select.pl -i 13,4 | sort -u | tee pbPA2TTS.tsv | cut -f1 | box.R -ng -nJ -no -y='Distance to TTS' -p=pbPA2TTS.pdf".format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True, executable="/bin/bash")

def drawSSmotif(asMotif=None, outFile=None):
    with open(asMotif) as f:
        lineList = f.readlines()
        mySum = sum([int(i.strip("\n").split("\t")[1]) for i in lineList])
        tmp = open("tmp.txt", "w")
        if len(lineList) < 4:
            for i in lineList:
                lineInfo = i.strip("\n").split("\t")
                print >> tmp, "\t".join([lineInfo[0], str(float(lineInfo[1])/mySum)])
        else:
            for i in lineList[0:3]:
                lineInfo = i.strip("\n").split("\t")
                print >> tmp, "\t".join([lineInfo[0], str(float(lineInfo[1]) / mySum)])
            otherSum = sum([int(i.strip("\n").split("\t")[1]) for i in lineList[3:]])
            print >> tmp, "\t".join(["Other", str(float(otherSum)/mySum)])
        tmp.close()
        cmd = "cat tmp.txt | bar.R -fillV=V1 -fp -lgPos=top -w=12 -p={}.ssMotif.pdf 2>/dev/null".format(outFile)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        os.remove("tmp.txt")
        os.remove(asMotif)

def splicesite2figure(refParams=None):
    refGenome = refParams.ref_genome
    cmd = '''splicesite2Seq.pl -f {} SE.inc.splicesite | cut -f5 | sort | uniq -c | tr -s ' ' | tr ' ' "\t" | sort -n -k1,1nr | select.pl --index 3,2 >SE.inc.tmp'''.format(refGenome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    drawSSmotif("SE.inc.tmp", "SE.inc")
    cmd = '''splicesite2Seq.pl -f {} SE.exc.splicesite | cut -f5 | sort | uniq -c | tr -s ' ' | tr ' ' "\t" | sort -n -k1,1nr | select.pl --index 3,2 >SE.exc.tmp'''.format(refGenome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    drawSSmotif("SE.exc.tmp", "SE.exc")
    cmd = '''splicesite2Seq.pl -f {} --no3 A5SS.inc.splicesite | cut -f5 | sort | uniq -c | tr -s ' ' | tr ' ' "\t" | sort -n -k1,1nr | select.pl --index 3,2 >A5SS.inc.tmp'''.format(refGenome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    drawSSmotif("A5SS.inc.tmp", "A5SS.inc")
    cmd = '''splicesite2Seq.pl -f {} --no3 A5SS.exc.splicesite | cut -f5 | sort | uniq -c | tr -s ' ' | tr ' ' "\t" | sort -n -k1,1nr | select.pl --index 3,2 >A5SS.exc.tmp'''.format(refGenome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    drawSSmotif("A5SS.exc.tmp", "A5SS.exc")
    cmd = '''splicesite2Seq.pl -f {} --no5 A3SS.inc.splicesite | cut -f5 | sort | uniq -c | tr -s ' ' | tr ' ' "\t" | sort -n -k1,1nr | select.pl --index 3,2 >A3SS.inc.tmp'''.format(refGenome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    drawSSmotif("A3SS.inc.tmp", "A3SS.inc")
    cmd = '''splicesite2Seq.pl -f {} --no5 A3SS.exc.splicesite | cut -f5 | sort | uniq -c | tr -s ' ' | tr ' ' "\t" | sort -n -k1,1nr | select.pl --index 3,2 >A3SS.exc.tmp'''.format(refGenome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    drawSSmotif("A3SS.exc.tmp", "A3SS.exc")
    cmd = '''splicesite2Seq.pl -f {} IR.splicesite | cut -f5 | sort | uniq -c | tr -s ' ' | tr ' ' "\t" | sort -n -k1,1nr | select.pl --index 3,2 >IR.tmp'''.format(refGenome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    drawSSmotif("IR.tmp", "IR")
    cmd = '''splicesite2Seq.pl -f {} IR.anno.splicesite | cut -f5 | sort | uniq -c | tr -s ' ' | tr ' ' "\t" | sort -n -k1,1nr | select.pl --index 3,2 >IR.anno.tmp'''.format(refGenome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    drawSSmotif("IR.anno.tmp", "IR.anno")
    cmd = '''splicesite2Seq.pl -f {} IR.novel.splicesite | cut -f5 | sort | uniq -c | tr -s ' ' | tr ' ' "\t" | sort -n -k1,1nr | select.pl --index 3,2 >IR.novel.tmp'''.format(refGenome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    drawSSmotif("IR.novel.tmp", "IR.novel")


def getSpliceSite(asType=None, asFile=None, outFile=None):
    if asType == "IR":
        out = open(outFile, "w")
        with open(asFile) as f:
            for line in f:
                lineInfo = line.strip("\n").split(":")
                posList = lineInfo[2].split("-")
                print >>out, "\t".join([lineInfo[0], posList[0], posList[1], lineInfo[3]])
        out.close()
    elif asType == "SE":
        cmd = "seDecompose.pl        confident.SE.lst   >SE.inc.splicesite   2>SE.exc.splicesite"
        subprocess.call(cmd, shell=True)
    elif asType == "A3SS":
        cmd = "anssDecompose.pl -n 5 confident.A5SS.lst >A5SS.inc.splicesite 2>A5SS.exc.splicesite"
        subprocess.call(cmd, shell=True)
    elif asType == "A5SS":
        cmd = "anssDecompose.pl -n 3 confident.A3SS.lst >A3SS.inc.splicesite 2>A3SS.exc.splicesite"
        subprocess.call(cmd, shell=True)

def identifyDE(researchObjs=None):
    # degDict = {"plotAll": 1, "genes": []}
    degDict = {}
    for i in researchObjs:
        paramObj = researchObjs[i]
        degDict[paramObj.project] = {"plotAll": 1, "genes": []}
        if paramObj.tgs_plat and paramObj.ngs_dir == None:
            if paramObj.tgs_plat.lower() == "pacbio":
                print "PacBio reads only can't be quantified! The structure of all the isoforms will be plotted!"
            elif paramObj.tgs_plat.lower() == "nanopore":
                print "Nanopore can be quantified! The differential expressed genes between samples will be plotted!"
        elif paramObj.tgs_plat and paramObj.ngs_dir:
            if paramObj.tgs_plat.lower() == "pacbio":
                print "It's PACBIO hybrid strategy! The differential expressed genes between samples will be plotted!"
            elif paramObj.tgs_plat.lower() == "nanopore":
                print "It's NANOPORE hybrid strategy! The differential expressed genes between samples will be plotted!"
        elif paramObj.tgs_plat == None and paramObj.ngs_dir:
            print "It's NGS strategy! The differential expressed genes between samples will be plotted!"
    return degDict

def isoViewer(refParams=None, sampleCfgFile=None, degDict=None):
    prevDir = os.getcwd()
    gpeObj = GenePredObj(refParams.ref_gpe, False)

    sampleDict = getSampleCfg(sampleCfgFile)
    gpeDict = gpeObj.genePredDict
    # geneName2gpeObj = gpeObj.geneName2gpeObj
    # sample2reads = {}
    sampleObj = Sample()
    for sample in sampleDict:
        sourceFile = sampleDict[sample]
        validateFile(sourceFile)
        with open(sourceFile) as f:
            for line in f:
                readStruc = ReadLineStruc(line)
                sampleObj.update(sample, sourceFile, readStruc)
                # getAS(gpeDict, readStruc, "", True, geneDict, novelDict, offset)
                # readStruc,gpeDict

        # testGene = "GRMZM2G336448"
        # testGene = "GRMZM2G056799"
        # testGene = "GRMZM2G029058"
        # testGene = "GRMZM2G008061"
        testGene = "GRMZM2G109009"
        geneDict, novelDict = {}, {}
        if degDict[sample]["plotAll"] == 0:
            pass
        else:
            pass
        for g in sampleObj.sample2Gene[sample]:
            getAS(gpeDict, sampleObj.sample2Gene[sample][g], geneDict, novelDict, offset=0)
            # if sampleObj.sample2Gene[sample][g].geneName == testGene:
            #     getAS(gpeDict, sampleObj.sample2Gene[sample][g], geneDict, novelDict, offset=0)
            # sample2reads[sample].update()
        geneDict = findAS(geneDict, anno=True)
        novelDict = findAS(novelDict, anno=False)

        targetDir = os.path.join(prevDir, sample, testGene)
        if os.path.exists(targetDir):
            shutil.rmtree(targetDir)
        os.makedirs(targetDir)
        os.chdir(targetDir)
        figOut = testGene + ".pdf"
        cfgOut = open(testGene + ".cfg", "w")
        mainSec = MainSection(fout=figOut)
        print >> cfgOut, mainSec.printStr()

        geneModelGPE = "%s.gpe" % (testGene)
        geneModelGTF = "%s.gtf" % (testGene)
        geneModelGFF = "%s.gff" % (testGene)
        gpeObj.gpeObj2file(geneModelGPE, testGene)
        os.system("genePredToGtf file %s %s" % (geneModelGPE, geneModelGTF))
        os.system("gene_model_to_splicegraph.py -m %s -o %s" % (geneModelGTF, geneModelGFF))
        geneModelHidePlot = PlotSection(section_name="[GeneModelGraph]", source_file=geneModelGTF, gene_name=testGene,
                                        title_string="Gene Model for %gene", hide=True)
        print >> cfgOut, geneModelHidePlot.printStr()
        geneModelVisiblePlot = PlotSection(section_name="[GeneModelIsoformsGraph]", plot_type="isoforms",
                                           source_file=geneModelGFF, title_string="Gene Model for %gene")
        print >> cfgOut, geneModelVisiblePlot.printStr()

        allReads = geneDict[testGene].reads
        allReadsGPE = "%s.allReads.gpe" % (testGene)
        allReadsGTF = "%s.allReads.gtf" % (testGene)
        allReadsGFF = "%s.allReads.gff" % (testGene)
        allReadsGPEOut = open(allReadsGPE, "w")
        for read in allReads:
            print >> allReadsGPEOut, allReads[read].go_to_gpe()
        allReadsGPEOut.close()
        os.system("genePredToGtf file %s %s" % (allReadsGPE, allReadsGTF))
        os.system("gene_model_to_splicegraph.py -m %s -o %s" % (allReadsGTF, allReadsGFF))
        allReadsPlot = PlotSection(section_name="[AllReadsCollapse]", plot_type="splice_graph",
                                   source_file=allReadsGFF, title_string="All Reads Presented in %s" % (testGene))
        print >> cfgOut, allReadsPlot.printStr()

        asDict = geneDict[testGene].asDict
        asType2asKey = {"A5SS": "exc", "A3SS": "exc", "SE": "skip", "IR": "spliced"}
        for asType in asDict:
            asGPE = "%s.%s.gpe" % (testGene, asType)
            asGTF = "%s.%s.gtf" % (testGene, asType)
            asGFF = "%s.%s.gff" % (testGene, asType)
            print asType
            if not asDict[asType]:
                continue
            asGPEOut = open(asGPE, "w")
            if asType not in asType2asKey:
                for read in asDict[asType]:
                    print >> asGPEOut, read.go_to_gpe()
            else:
                excValue = asType2asKey[asType]
                for junc in asDict[asType]:
                    for read in asDict[asType][junc][excValue]:
                        print >> asGPEOut, read.go_to_gpe()
            asGPEOut.close()
            os.system("genePredToGtf file %s %s" % (asGPE, asGTF))
            os.system("gene_model_to_splicegraph.py -m %s -o %s -a" % (asGTF, asGFF))
            asPlot = PlotSection(section_name="[%s]" % (asType), plot_type="splice_graph", source_file=asGFF,
                                 title_string="%s in %s" % (asType, testGene))
            print >> cfgOut, asPlot.printStr()

        # coverage
        # readsPlot = PlotSection(section_name="[Reads]", plot_type="read_depth", source_file="../../../data/uniq.sorted.sam",
        #                         title_string="%s Read Coverage" % (testGene))
        # print >>cfgOut, readsPlot.printStr()
        cfgOut.close()
        os.system("plotter.py %s" % (testGene + ".cfg"))
        os.chdir(prevDir)

def singleRun(refParams=None, paramObj=None, ccsParams=None, mecatParams=None, project2as=None):
    initSysSetting(refParams=refParams)
    projectName = paramObj.project
    prevDir = os.getcwd()
    workDir = os.path.join(refParams.out_dir, projectName)
    if os.path.exists(workDir):
        shutil.rmtree(workDir)
    os.mkdir(workDir)
    os.chdir(workDir)
    processRawReads(paramObj=paramObj, refParams=refParams, ccsParams=ccsParams, mecatParams=mecatParams)
    # makeRefTranscriptome()
    filterAndRefineAlign(refParams=refParams)
    removeArtifacts(refParams=refParams)
    collapse(tgs_plat=paramObj.tgs_plat, refParams=refParams)
    postMappingEvaluation(refParams=refParams)
    paAnalysis(refParams=refParams)
    identifyASE(refParams=refParams, paramObj=paramObj, project2as=project2as)
    os.chdir(prevDir)

def multiTest(refParams=None, paramObj=None, ccsParams=None, mecatParams=None, project2as=None):
    b = []
    for i in range(20000000):
        b.append(i)
    return "multiTest done! {}, {}, {}, {}, {}".format(refParams, paramObj, ccsParams, mecatParams, project2as)

def multiRun(defaultCfg=None):
    researchObjs = defaultCfg.researchObjs
    refParams = defaultCfg.refParams
    mecatParams = defaultCfg.mecatParams
    ccsParams = defaultCfg.ccsParams
    project2as = {}
    # project2as = {"project1": "/home/xufeng/xufeng/iso-seq/test/workflow_scripts_test/test_20190808/test/project1/ASE/deSingleExonRead.bed12+",
    #               "project2": "/home/xufeng/xufeng/iso-seq/test/workflow_scripts_test/test_20190808/test/project2/ASE/deSingleExonRead.bed12+"}
    researchNum = len(researchObjs)
    pool = Pool(processes=researchNum)
    multiResults = []
    for i in researchObjs:
        singleRunRes = pool.apply_async(singleRun, (refParams, researchObjs[i], ccsParams, mecatParams, project2as))
        multiResults.append(singleRunRes)
    for j in multiResults:
        j.wait()
    for j in multiResults:
        if j.ready():
            if j.successful():
                print j.get()
    # for i in researchObjs:
    #     singleRun(refParams=refParams, paramObj=researchObjs[i], ccsParams=ccsParams, mecatParams=mecatParams, project2as=project2as)

    degDict = identifyDE(researchObjs=researchObjs)
    sampleOut = open("sample.conf", "w")
    for i in project2as:
        print >>sampleOut, "\t".join([i, project2as[i]])
    sampleOut.close()
    isoViewer(refParams=refParams, sampleCfgFile="sample.conf", degDict=degDict)

def main():
    if args.command:
        paramObj = ResearchSection("singleRun")
        defaultCfg = None
        for key, value in args._get_kwargs():
            if key == "default_cfg":
                defaultCfg = Config(args.default_cfg)
            else:
                paramObj.__dict__[key] = value
        print paramObj, defaultCfg

        refParams = defaultCfg.refParams
        mecatParams = defaultCfg.mecatParams
        ccsParams = defaultCfg.ccsParams
        singleRun(refParams=refParams, paramObj=paramObj, ccsParams=ccsParams, mecatParams=mecatParams)
    else:
        defaultCfg = Config(args.default_cfg)
        multiRun(defaultCfg=defaultCfg)

if __name__ == '__main__':
    main()
