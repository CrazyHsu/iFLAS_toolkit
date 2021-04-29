#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: identifyASE.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-05 15:55:41
Last modified: 2019-12-05 15:55:41
'''

import datetime, argparse, pybedtools
from multiprocessing import Pool

from commonFuncs import *
from commonObjs import *
# from findAS import *
from Config import *
from identifyASEFuncs import *
from scanESbyNGSFuncs import *
from assignGeneNameToSE import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

##############################################
def identifyASE(refParams=None, tgsSample=None, dirSpec=None):
    # Identify alternative splicing events
    print str(datetime.datetime.now()) + " Alternative splicing events identifying for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    groupDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName)
    aseDir = os.path.join(groupDir, "ASE")
    resolveDir(aseDir)
    makeLink(os.path.join(groupDir, "postCorrEval", "isoformGrouped.bed12+"), "isoformGrouped.bed12+")
    makeLink(os.path.join(groupDir, "filtration", "collapse", "tofu.collapsed.group.txt"), "tofu.collapsed.group.txt")
    if tgsSample.ngsJunctions == None:
        tgsSample.ngsJunctions = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "RNA-seq", "reassembly", "junctions.bed12")

    if not os.path.exists("PB"):
        os.makedirs("PB")
    transBedList = GenePredObj(refParams.ref_gpe, bincolumn=False).toBed(gene=True)
    transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)
    readsBedList = []
    isoform2reads = {}
    with open("tofu.collapsed.group.txt") as f:
        for line in f.readlines():
            isoform, reads = line.strip("\n").split("\t")
            isoform2reads[isoform] = reads.split(",")
    with open("isoformGrouped.bed12+") as f:
        for line in f.readlines():
            readStruc = Bed12(line)
            if len(isoform2reads[readStruc.name]) < 2: continue
            readsBedList.append("\t".join(readStruc.record[:12]))
    readsBedObj = pybedtools.BedTool("\n".join(readsBedList), from_string=True)

    annoBedRes = readsBedObj.intersect(transBedObj, wa=True, wb=True, s=True)
    novelBedRes = readsBedObj.intersect(transBedObj, v=True, s=True)
    annoDict, novelDict = getAS(annoBedRes, novelBedRes, offset=0)
    findASmain_bak(asType="IR", annoDict=annoDict, novelDict=novelDict, outFile="PB/IR.bed6+", isoform2reads=isoform2reads)
    findASmain_bak(asType="SE", annoDict=annoDict, novelDict=novelDict, outFile="PB/SE.bed12+", isoform2reads=isoform2reads)
    findASmain_bak(asType="A5SS", annoDict=annoDict, novelDict=novelDict, outFile="PB/A5SS.bed6+", isoform2reads=isoform2reads)
    findASmain_bak(asType="A3SS", annoDict=annoDict, novelDict=novelDict, outFile="PB/A3SS.bed6+", isoform2reads=isoform2reads)

    # findASmain(asType="IR", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/IR.bed6+")
    # findASmain(asType="SE", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/SE.bed12+")
    # findASmain(asType="A5SS", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/A5SS.bed6+")
    # findASmain(asType="A3SS", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/A3SS.bed6+")


    # cmd = "getSingleIr.pl -g {} deSingleExonRead.bed12+ >PB/IR.bed6+".format(refParams.ref_gpe)
    # subprocess.call(cmd, shell=True)
    # cmd = "getSE.pl -g {} deSingleExonRead.bed12+ >PB/SE.bed12+".format(refParams.ref_gpe)
    # subprocess.call(cmd, shell=True)
    # cmd = "getA5SS.pl -g {} deSingleExonRead.bed12+ >PB/A5SS.bed6+".format(refParams.ref_gpe)
    # subprocess.call(cmd, shell=True)
    # cmd = "getA3SS.pl -g {} deSingleExonRead.bed12+ >PB/A3SS.bed6+".format(refParams.ref_gpe)
    # subprocess.call(cmd, shell=True)

    if os.path.islink("PB/IR.confident.bed6+") or os.path.exists("PB/IR.confident.bed6+"):
        os.remove("PB/IR.confident.bed6+")
    if os.path.islink("PB/SE.confident.bed12+") or os.path.exists("PB/SE.confident.bed12+"):
        os.remove("PB/SE.confident.bed12+")
    if os.path.islink("PB/A5SS.confident.bed12+") or os.path.exists("PB/A5SS.confident.bed6+"):
        os.remove("PB/A5SS.confident.bed6+")
    if os.path.islink("PB/A3SS.confident.bed6+") or os.path.exists("PB/A3SS.confident.bed6+"):
        os.remove("PB/A3SS.confident.bed6+")

    if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
        if not os.path.exists("NGS"):
            os.makedirs("NGS")
        try:
            from scanESbyNGS import scanEsByNGS
            scanEsByNGS(fref=refParams.ref_gpe, has_bin=False, fjunc=tgsSample.ngsJunctions, fpb="isoformGrouped.bed12+", outFile="SEsupportByNGS.bed12+")
            assignGeneNameToSE(fref=refParams.ref_gpe, has_bin=False, fes="SEsupportByNGS.bed12+", outFile="NGS/SE.bed12+", errorOutFile="NGS/novelSE2gene.log")
        except Exception as e:
            print e
        filter(originFile="NGS/SE.bed12+", targetFile="PB/SE.bed12+", originField=4, targetField=4,
               outFile="PB/SE.NGS.bed12+")
        filter(originFile="PB/SE.bed12+", targetFile="NGS/SE.bed12+", originField=4, targetField=4,
               outFile="NGS/SE.PB.bed12+")
        cmd = "juncAssign.pl -g {} {} | tee junction.assigned.for_IR.bed12+ | cut -f 1-12,21-22 >junction.assigned.bed12+".format(refParams.ref_gpe, tgsSample.ngsJunctions)
        subprocess.call(cmd, shell=True)
        cmd = '''
                AnSSconfirmByJunc.pl -t 5 -j junction.assigned.bed12+ PB/A5SS.bed6+ >NGS/A5SS.PB.bed6+ 2>PB/A5SS.NGS.bed6+
                AnSSconfirmByJunc.pl -t 3 -j junction.assigned.bed12+ PB/A3SS.bed6+ >NGS/A3SS.PB.bed6+ 2>PB/A3SS.NGS.bed6+
        '''
        subprocess.call(cmd, shell=True)

        # irConfirmByJunc("PB/IR.bed6+", "junction.assigned.bed12+", ngsOutFile="NGS/IR.PB.bed6+", pbOutFile="PB/IR.NGS.bed6+")
#        irConfirmByJunc_1("PB/IR.bed6+", "junction.assigned.for_IR.bed12+", offset=3, juncBedAssigned=True,
#                          ngsOutFile="NGS/IR.PB.bed6+", pbOutFile="PB/IR.NGS.bed6+")

        cmd = '''
                awk '$9>1 && $11>1{print $3-$2}' NGS/A5SS.PB.bed6+ | hist.R -p=NGS/A5SS.PB.size.pdf 2>/dev/null
                awk '$9>1 && $11>1{print $5/1000}' NGS/A5SS.PB.bed6+ | hist.R -p=NGS/A5SS.PB.InclusionRatio.pdf 2>/dev/null
                awk '$9>1 && $11>1{print $3-$2}' NGS/A3SS.PB.bed6+ | hist.R -p=NGS/A3SS.PB.size.pdf 2>/dev/null
                awk '$9>1 && $11>1{print $5/1000}' NGS/A3SS.PB.bed6+ | hist.R -p=NGS/A3SS.PB.InclusionRatio.pdf 2>/dev/null
        '''
        subprocess.call(cmd, shell=True)
        cmd = '''
                awk 'BEGIN{FS=OFS="\t"}{$4=$4":"$6;print}' junction.assigned.bed12+ | getNgsA5SS.pl -e 0 >NGS/A5SS.known.bed6+ 2>NGS/A5SS.novel.bed6+
                awk 'BEGIN{FS=OFS="\t"}{$4=$4":"$6;print}' junction.assigned.bed12+ | getNgsA3SS.pl -e 0 >NGS/A3SS.known.bed6+ 2>NGS/A3SS.novel.bed6+
        '''
        ## need add IR and SE NGS annotation
        subprocess.call(cmd, shell=True)
        cmd = '''
                awk '$9>1 && $11>1{print $3-$2}' NGS/A5SS.{known,novel}.bed6+ | hist.R -p=NGS/A5SS.size.pdf 2>/dev/null
                awk '$9>1 && $11>1{print $5/1000}' NGS/A5SS.{known,novel}.bed6+ | hist.R -p=NGS/A5SS.InclusionRatio.pdf 2>/dev/null
                awk '$9>1 && $11>1{print $3-$2}' NGS/A3SS.{known,novel}.bed6+ | hist.R -p=NGS/A3SS.size.pdf 2>/dev/null
                awk '$9>1 && $11>1{print $5/1000}' NGS/A3SS.{known,novel}.bed6+ | hist.R -p=NGS/A3SS.InclusionRatio.pdf 2>/dev/null
        '''
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        makeLink("IR.NGS.bed6+", "PB/IR.confident.bed6+")
        makeLink("SE.NGS.bed12+", "PB/SE.confident.bed12+")
        makeLink("A5SS.NGS.bed6+", "PB/A5SS.confident.bed6+")
        makeLink("A3SS.NGS.bed6+", "PB/A3SS.confident.bed6+")
    else:
        makeLink("IR.bed+", "PB/IR.confident.bed12+")
        makeLink("SE.bed12+", "PB/SE.confident.bed12+")
        makeLink("A5SS.bed6+", "PB/A5SS.confident.bed6+")
        makeLink("A3SS.bed6+", "PB/A3SS.confident.bed6+")

    cmd = '''(cut -f 8,10 --output-delimiter=',' PB/A3SS.confident.bed6+ PB/A5SS.confident.bed6+ PB/IR.confident.bed6+;
              cut -f 16,18 --output-delimiter=',' PB/SE.confident.bed12+) | tr ',' '\n' | sort -u |
              filter.pl -o - isoformGrouped.bed12+ -2 4 -m i > isoformGrouped.AS.confident.bed12+'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    # if pbSample.projectName not in project2as:
    #     project2as[pbSample.projectName] = {pbSample.group: os.path.join(os.getcwd(), "deSingleExonIsoform.AS.confident.bed12+")}
    # else:
    #     project2as[pbSample.projectName].update({pbSample.group: os.path.join(os.getcwd(), "deSingleExonIsoform.AS.confident.bed12+")})

    cmd = '''
            awk '{print $3-$2}' PB/A5SS.confident.bed6+ | hist.R -p=PB/A5SS.size.pdf 2>/dev/null
            awk '$9>1 && $11>1{print $5/1000}' PB/A5SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=PB/A5SS.InclusionRatio.pdf 2>/dev/null
            awk '{print $3-$2}' PB/A3SS.confident.bed6+ | hist.R -p=PB/A3SS.size.pdf 2>/dev/null
            awk '$9>1 && $11>1{print $5/1000}' PB/A3SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=PB/A3SS.InclusionRatio.pdf 2>/dev/null
        '''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Alternative splicing events identifying for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)
    characterizeASE(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)

def characterizeASE(refParams=None, tgsSample=None, dirSpec=None):
    # Alternative splicing events characterization
    print str(datetime.datetime.now()) + " ASE Characterization for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    characAseDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "ASE", "characterization")
    resolveDir(characAseDir)
    getAnnoASList("../PB/SE.bed12+", "PB.SE.lst")
    getAnnoASList("../PB/A5SS.bed6+", "PB.A5SS.lst")
    getAnnoASList("../PB/A3SS.bed6+", "PB.A3SS.lst")
    getAnnoASList("../PB/IR.bed6+", "PB.IR.lst")
    getAnnoASList("../../PA/reads.paGrouped.tsv", "PB.APA.lst", PA=True)

    if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
        getAnnoASList("../NGS/SE.bed12+", "NGS.SE.lst")
        getAnnoASList("../NGS/A5SS.known.bed6+", "NGS.A5SS.lst")
        getAnnoASList("../NGS/A5SS.PB.bed6+", "NGS.A5SS.lst", append=True, uniq=True)
        getAnnoASList("../NGS/A3SS.known.bed6+", "NGS.A3SS.lst")
        getAnnoASList("../NGS/A3SS.PB.bed6+", "NGS.A3SS.lst", append=True, uniq=True)

        filter(originFile="NGS.SE.lst", targetFile="PB.SE.lst", outFile="Common.SE.lst")
        filter(originFile="NGS.A5SS.lst", targetFile="PB.A5SS.lst", outFile="Common.A5SS.lst")
        filter(originFile="NGS.A3SS.lst", targetFile="PB.A3SS.lst", outFile="Common.A3SS.lst")

        makeLink("Common.SE.lst", "confident.SE.lst")
        makeLink("Common.A5SS.lst", "confident.A5SS.lst")
        makeLink("Common.A3SS.lst", "confident.A3SS.lst")

        cmd = '''
                venn.R {NGS,PB}.SE.lst   -rD=90 -cP=180,0 -p=SE.pdf 1>/dev/null 2>&1
                venn.R {NGS,PB}.A5SS.lst -rD=90 -cP=180,0 -p=A5SS.pdf 1>/dev/null 2>&1
                venn.R {NGS,PB}.A3SS.lst -rD=90 -cP=180,0 -p=A3SS.pdf 1>/dev/null 2>&1
        '''
        subprocess.call(cmd, shell=True)

        filter(originFile="NGS.SE.lst", targetFile="PB.SE.lst", outFile="PB.specific.SE.lst")
        filter(originFile="NGS.A5SS.lst", targetFile="PB.A5SS.lst", outFile="PB.specific.A5SS.lst")
        filter(originFile="NGS.A3SS.lst", targetFile="PB.A3SS.lst", outFile="PB.specific.A3SS.lst")

        cmd = "seDiagnose.pl -j {} PB.specific.SE.lst > seDiagnose.tsv".format(tgsSample.ngsJunctions)
        subprocess.call(cmd, shell=True)
        cmd = "anssDiagnose.pl -a 5 -j {} PB.specific.A5SS.lst >a5ssDiagnose.tsv".format(tgsSample.ngsJunctions)
        subprocess.call(cmd, shell=True)
        cmd = "anssDiagnose.pl -a 3 -j {} PB.specific.A3SS.lst >a3ssDiagnose.tsv".format(tgsSample.ngsJunctions)
        subprocess.call(cmd, shell=True)
    else:
        makeLink("PB.SE.lst", "confident.SE.lst")
        makeLink("PB.A5SS.lst", "confident.A5SS.lst")
        makeLink("PB.A3SS.lst", "confident.A3SS.lst")

    # Classify APE into annotated or novel according to reference gene model
    cmd = '''awk '{{print $0"\t"$12}}' {} | gpe2bed.pl -p >reference.bed12+'''.format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    transBedList = GenePredObj(refParams.ref_gpe, bincolumn=False).toBed(gene=True)
    transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)
    readsBedList = []
    with open("reference.bed12+") as f:
        for line in f.readlines():
            readStruc = Bed12(line)
            readsBedList.append("\t".join(readStruc.record[:12]))
    readsBedObj = pybedtools.BedTool("\n".join(readsBedList), from_string=True)

    annoBedRes = readsBedObj.intersect(transBedObj, wa=True, wb=True, s=True)
    novelBedRes = readsBedObj.intersect(transBedObj, v=True, s=True)
    annoDict, novelDict = getAS(annoBedRes, novelBedRes, offset=0)
    findASmain_bak(asType="IR", annoDict=annoDict, novelDict=novelDict, outFile="IR.reference.bed6+")
    findASmain_bak(asType="SE", annoDict=annoDict, novelDict=novelDict, outFile="SE.reference.bed6+")
    findASmain_bak(asType="A5SS", annoDict=annoDict, novelDict=novelDict, outFile="A5SS.reference.bed6+")
    findASmain_bak(asType="A3SS", annoDict=annoDict, novelDict=novelDict, outFile="A3SS.reference.bed6+")

    # findASmain(asType="IR", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="IR.reference.bed6+")
    # findASmain(asType="SE", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="SE.reference.bed6+")
    # findASmain(asType="A3SS", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="A3SS.reference.bed6+")
    # findASmain(asType="A5SS", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="A5SS.reference.bed6+")

    # cmd = "getSingleIr.pl -g {} reference.bed12+ >IR.reference.bed6+".format(refParams.ref_gpe)
    # subprocess.call(cmd, shell=True)
    # cmd = "getSE.pl -g {} reference.bed12+ >SE.reference.bed6+".format(refParams.ref_gpe)
    # subprocess.call(cmd, shell=True)
    # cmd = "getA5SS.pl -g {} reference.bed12+ >A5SS.reference.bed6+".format(refParams.ref_gpe)
    # subprocess.call(cmd, shell=True)
    # cmd = "getA3SS.pl -g {} reference.bed12+ >A3SS.reference.bed6+".format(refParams.ref_gpe)
    # subprocess.call(cmd, shell=True)
    getAnnoASList("IR.reference.bed6+", "IR.reference.lst")
    getAnnoASList("SE.reference.bed6+", "SE.reference.lst")
    getAnnoASList("A3SS.reference.bed6+", "A3SS.reference.lst")
    getAnnoASList("A5SS.reference.bed6+", "A5SS.reference.lst")
    cmd = "paGroup.pl reference.bed12+ >paGroup.reference.tsv 2>paGroup.reference.bed6"
    subprocess.call(cmd, shell=True)
    cmd = "paCmp.pl -r paGroup.reference.tsv -a annoPA.tsv -n novelPA.tsv ../../PA/reads.paGrouped.tsv >paGroup.novel.tsv 2>paGroup.anno.tsv"
    subprocess.call(cmd, shell=True)

    filter(originFile="IR.reference.lst", targetFile="PB.IR.lst", outFile="IR.anno.lst", mode="i")
    filter(originFile="IR.reference.lst", targetFile="PB.IR.lst", outFile="IR.novel.lst", mode="e")
    filter(originFile="SE.reference.lst", targetFile="confident.SE.lst", outFile="SE.anno.lst", mode="i")
    filter(originFile="SE.reference.lst", targetFile="confident.SE.lst", outFile="SE.novel.lst", mode="e")
    filter(originFile="A3SS.reference.lst", targetFile="confident.A3SS.lst", outFile="A3SS.anno.lst", mode="i")
    filter(originFile="A3SS.reference.lst", targetFile="confident.A3SS.lst", outFile="A3SS.novel.lst", mode="e")
    filter(originFile="A5SS.reference.lst", targetFile="confident.A5SS.lst", outFile="A5SS.anno.lst", mode="i")
    filter(originFile="A5SS.reference.lst", targetFile="confident.A5SS.lst", outFile="A5SS.novel.lst", mode="e")

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

    # Event Decompose and Motif Evaluation
    getSpliceSite(asType="IR", asFile="PB.IR.lst", outFile="IR.splicesite")
    getSpliceSite(asType="IR", asFile="IR.anno.lst", outFile="IR.anno.splicesite")
    getSpliceSite(asType="IR", asFile="IR.novel.lst", outFile="IR.novel.splicesite")
    getSpliceSite(asType="SE", asFile="confident.SE.lst")
    getSpliceSite(asType="A3SS", asFile="confident.A3SS.lst")
    getSpliceSite(asType="A5SS", asFile="confident.A5SS.lst")
    splicesite2figure(refParams=refParams)

    getDist2TTS(refParams=refParams, paGroup="../../PA/reads.paGrouped.bed6")

    if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
        cmd = '''
            skyjoin <(cut -f4,5 ../PB/SE.bed12+) <(cut -f4,5 ../NGS/SE.bed12+) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=1-$3/1000;print $2,$3}' \
            | correlation.R -x=PB -y=NGS -p=SE.cor.pdf >SE.cor 2>/dev/null
            skyjoin <(cut -f4,5 ../PB/SE.bed12+) <(cut -f4,5 ../NGS/SE.bed12+) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=1-$3/1000;print "PB",$2,$1;print "NGS",$3,$1}' \
            | line.R -b=V3 -s=0.01 -a=0.05 -x=Platform -y=PSI -m="R:$(head -n1 SE.cor | cut -f2), P:$(head -n1 SE.cor | cut -f3)" -p=SE.cor2.pdf 2>/dev/null
            skyjoin <(cut -f4,5 ../PB/A5SS.bed6+) <(cut -f4,5 ../NGS/A5SS.{known,PB}.bed6+ | sort -u) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print $2,$3}' \
            | correlation.R -x=PB -y=NGS -p=A5SS.cor.pdf >A5SS.cor 2>/dev/null
            skyjoin <(cut -f4,5 ../PB/A5SS.bed6+) <(cut -f4,5 ../NGS/A5SS.{known,PB}.bed6+ | sort -u) 2>/dev/null \
            | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print "PB",$2,$1;print "NGS",$3,$1}' \
            | line.R -b=V3 -s=0.01 -a=0.05 -x=Platform -y=PSI -m="R:$(head -n1 A5SS.cor | cut -f2), P:$(head -n1 A5SS.cor | cut -f3)" -p=A5SS.cor2.pdf 2>/dev/null
            skyjoin <(cut -f4,5 ../PB/A3SS.bed6+) <(cut -f4,5 ../NGS/A3SS.{known,PB}.bed6+ | sort -u) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print $2,$3}' \
            | correlation.R -x=PB -y=NGS -p=A3SS.cor.pdf >A3SS.cor 2>/dev/null
            skyjoin <(cut -f4,5 ../PB/A3SS.bed6+) <(cut -f4,5 ../NGS/A3SS.{known,PB}.bed6+ | sort -u) 2>/dev/null \
            | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print "PB",$2,$1;print "NGS",$3,$1}' \
            | line.R -b=V3 -s=0.01 -a=0.05 -x=Platform -y=PSI -m="R:$(head -n1 A3SS.cor | cut -f2), P:$(head -n1 A3SS.cor | cut -f3)" -p=A3SS.cor2.pdf 2>/dev/null
        '''
        subprocess.call(cmd, shell=True, executable="/bin/bash")
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " ASE Characterization for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)

# def identifyASE_bak(refParams=None, pbSample=None, project2as=None):
#     # Identify alternative splicing events
#     print str(datetime.datetime.now()) + " Alternative splicing events identifying for group {}...".format(pbSample.group)
#     aseDir = os.path.join(refParams.out_dir, pbSample.projectName, pbSample.group, "ASE")
#     resolveDir(aseDir)
#     cmd = "readGroup.pl -g {} ../postCorrEval/assigned.bed12+ >unambi.bed12+ 2>ambiguous.bed12+".format(refParams.ref_gpe)
#     subprocess.call(cmd, shell=True)
#     with open("unambi.bed12+") as f:
#         deSingleExonRead = open("deSingleExonRead.bed12+", "w")
#         readGroup = open("readGrouped.bed12+", "w")
#         meDict = getColFromFile("../PA/singleExonReadWithExonInMEread.bed12+", 4, sep=",")
#         seDict = getColFromFile("../PA/paCluster.singleExonRead.bed8+", 4, sep=",")
#         for line in f:
#             lineInfo = line.strip("\n").split("\t")
#             if lineInfo[3] not in meDict and lineInfo[3] not in seDict:
#                 print >>deSingleExonRead, "\t".join(lineInfo[0:12]+lineInfo[20:21])
#             print >>readGroup, "\t".join(lineInfo[0:12]+lineInfo[20:21])
#         deSingleExonRead.close()
#         readGroup.close()
#     # cmd = '''cut -f1-12,21 unambi.bed12+ | tee readGrouped.bed12+ | filter.pl -o <(cut -f4 ../PA/singleExonReadWithExonInMEread.bed12+; cut -f4 ../PA/paCluster.singleExonRead.bed8+ | tr ',' "\n") -2 4 >deSingleExonRead.bed12+'''
#     cmd = "3endRevise.pl -p ../PA/PA.bed6+ deSingleExonRead.bed12+ | tee deSingleExonRead.3endRevised.bed12+ | paGroup.pl >paGrouped.tsv 2>paGrouped.bed6"
#     subprocess.call(cmd, shell=True)
#     if not os.path.exists("PB"):
#         os.makedirs("PB")
#     findASmain(asType="IR", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/IR.bed6+")
#     # cmd = "getSingleIr.pl -g {} deSingleExonRead.bed12+ >PB/IR.bed6+".format(refParams.ref_gpe)
#     # subprocess.call(cmd, shell=True)
#     geneWithIR = len(getColFromFile("PB/IR.bed6+", 7))
#     deSingleExonReadCount = len(getColFromFile("deSingleExonRead.bed12+", 13))
#     geneWithIrRatio = round(float(geneWithIR)/deSingleExonReadCount, 4)
#     findASmain(asType="SE", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/SE.bed12+")
#     findASmain(asType="A5SS", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/A5SS.bed6+")
#     findASmain(asType="A3SS", gpeFile=refParams.ref_gpe, sourceFile="deSingleExonRead.bed12+", outFile="PB/A3SS.bed6+")
#
#     if pbSample.projectName not in project2as:
#         project2as[pbSample.projectName] = {pbSample.group: os.path.join(os.getcwd(), "deSingleExonRead.bed12+")}
#     else:
#         project2as[pbSample.projectName].update({pbSample.group: os.path.join(os.getcwd(), "deSingleExonRead.bed12+")})
#     # cmd = "getSE.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/SE.bed12+ && getA5SS.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/A5SS.bed6+ && getA3SS.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/A3SS.bed6+"
#
#     cmd = '''awk '$9>1 && $11>1 && $5>=100 && $5<=900{print $3-$2}' PB/A5SS.bed6+ | hist.R -p=PB/A5SS.size.pdf 2>/dev/null'''
#     subprocess.call(cmd, shell=True)
#     cmd = '''awk '$3-$2>=20 && $9>1 && $11>1{print $5/1000}' PB/A5SS.bed6+ | hist.R -x='Incusion Ratio' -p=PB/A5SS.InclusionRatio.pdf 2>/dev/null'''
#     subprocess.call(cmd, shell=True)
#
#     if os.path.islink("PB/SE.confident.bed12+") or os.path.exists("PB/SE.confident.bed12+"):
#         os.remove("PB/SE.confident.bed12+")
#     if os.path.islink("PB/SE.confident.bed12+") or os.path.exists("PB/A5SS.confident.bed6+"):
#         os.remove("PB/A5SS.confident.bed6+")
#     if os.path.islink("PB/A3SS.confident.bed6+") or os.path.exists("PB/A3SS.confident.bed6+"):
#         os.remove("PB/A3SS.confident.bed6+")
#
#     if pbSample.ngsProject != None:
#         if not os.path.exists("NGS"):
#             os.makedirs("NGS")
#         from scanESbyNGS import scanEsByNGS
#         scanEsByNGS(fref=refParams.ref_gpe, has_bin=False, fjunc=pbSample.ngsJunctions, fpb="deSingleExonRead.bed12+", outFile="SEsupportByNGS.bed12+")
#         assignGeneNameToSE(fref=refParams.ref_gpe, has_bin=False, fes="SEsupportByNGS.bed12+", outFile="NGS/SE.bed12+")
#         filter(originFile="NGS/SE.bed12+", targetFile="PB/SE.bed12+", originField=4, targetField=4,
#                outFile="PB/SE.NGS.bed12+")
#         filter(originFile="PB/SE.bed12+", targetFile="NGS/SE.bed12+", originField=4, targetField=4,
#                outFile="NGS/SE.PB.bed12+")
#         cmd = "juncAssign.pl -g {} {} >junction.assigned.bed12+".format(refParams.ref_gpe, pbSample.ngsJunctions)
#         subprocess.call(cmd, shell=True)
#         cmd = '''
#                 AnSSconfirmByJunc.pl -t 5 -j junction.assigned.bed12+ PB/A5SS.bed6+ >NGS/A5SS.PB.bed6+ 2>PB/A5SS.NGS.bed6+
#                 AnSSconfirmByJunc.pl -t 3 -j junction.assigned.bed12+ PB/A3SS.bed6+ >NGS/A3SS.PB.bed6+ 2>PB/A3SS.NGS.bed6+
#         '''
#         subprocess.call(cmd, shell=True)
#         cmd = '''
#                 awk '$9>1 && $11>1{print $3-$2}' NGS/A5SS.PB.bed6+ | hist.R -p=NGS/A5SS.PB.size.pdf 2>/dev/null
#                 awk '$9>1 && $11>1{print $5/1000}' NGS/A5SS.PB.bed6+ | hist.R -p=NGS/A5SS.PB.InclusionRatio.pdf 2>/dev/null
#                 awk 'BEGIN{FS=OFS="\t"}{$4=$4":"$6;print}' junction.assigned.bed12+ | getNgsA5SS.pl -e 0 >NGS/A5SS.known.bed6+ 2>NGS/A5SS.novel.bed6+
#                 awk 'BEGIN{FS=OFS="\t"}{$4=$4":"$6;print}' junction.assigned.bed12+ | getNgsA3SS.pl -e 0 >NGS/A3SS.known.bed6+ 2>NGS/A3SS.novel.bed6+
#                 awk '$9>1 && $11>1{print $3-$2}' NGS/A5SS.{known,novel}.bed6+ | hist.R -p=NGS/A5SS.size.pdf 2>/dev/null
#                 awk '$9>1 && $11>1{print $5/1000}' NGS/A5SS.{known,novel}.bed6+ | hist.R -p=NGS/A5SS.InclusionRatio.pdf 2>/dev/null
#         '''
#         subprocess.call(cmd, shell=True)
#         makeLink("SE.NGS.bed12+", "PB/SE.confident.bed12+")
#         makeLink("A5SS.NGS.bed6+", "PB/A5SS.confident.bed6+")
#         makeLink("A3SS.NGS.bed6+", "PB/A3SS.confident.bed6+")
#     else:
#         makeLink("SE.bed12+", "PB/SE.confident.bed12+")
#         makeLink("A5SS.bed6+", "PB/A5SS.confident.bed6+")
#         makeLink("A3SS.bed6+", "PB/A3SS.confident.bed6+")
#
#     # cmd = "ln -sf SE.bed12+ PB/SE.confident.bed12+ && ln -sf A5SS.bed6+ PB/A5SS.confident.bed6+ && ln -sf A3SS.bed6+ PB/A3SS.confident.bed6+"
#     cmd = '''
#             awk '{print $3-$2}' PB/A5SS.confident.bed6+ | hist.R -p=PB/A5SS.size.pdf 2>/dev/null
#             awk '$9>1 && $11>1{print $5/1000}' PB/A5SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=PB/A5SS.InclusionRatio.pdf 2>/dev/null
#             awk '{print $3-$2}' PB/A3SS.confident.bed6+ | hist.R -p=PB/A3SS.size.pdf 2>/dev/null
#             awk '$9>1 && $11>1{print $5/1000}' PB/A3SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=PB/A3SS.InclusionRatio.pdf 2>/dev/null
#         '''
#     subprocess.call(cmd, shell=True, executable="/bin/bash")
#     print str(datetime.datetime.now()) + " Alternative splicing events identifying for group {} done!".format(pbSample.group)
#     characterizeASE_bak(refParams=refParams, pbSample=pbSample)
#
# def characterizeASE_bak(refParams=None, pbSample=None):
#     # Alternative splicing events characterization
#     print str(datetime.datetime.now()) + " ASE Characterization for group {}...".format(pbSample.group)
#     prevDir = os.getcwd()
#     resolveDir("characterization")
#     getAnnoASList("../PB/SE.bed12+", "PB.SE.lst")
#     getAnnoASList("../PB/A5SS.bed6+", "PB.A5SS.lst")
#     getAnnoASList("../PB/A3SS.bed6+", "PB.A3SS.lst")
#     getAnnoASList("../PB/IR.bed6+", "PB.IR.lst")
#     # getAnnoASList("../paGrouped.tsv", "PB.APA.lst", PA=True)
#
#     if pbSample.ngsProject != None:
#         getAnnoASList("../NGS/SE.bed12+", "NGS.SE.lst")
#         getAnnoASList("../NGS/A5SS.known.bed6+", "NGS.A5SS.lst")
#         getAnnoASList("../NGS/A5SS.PB.bed6+", "NGS.A5SS.lst", append=True, uniq=True)
#         getAnnoASList("../NGS/A3SS.known.bed6+", "NGS.A3SS.lst")
#         getAnnoASList("../NGS/A5SS.PB.bed6+", "NGS.A3SS.lst", append=True, uniq=True)
#
#         filter(originFile="NGS.SE.lst", targetFile="PB.SE.lst", outFile="Common.SE.lst")
#         filter(originFile="NGS.A5SS.lst", targetFile="PB.A5SS.lst", outFile="Common.A5SS.lst")
#         filter(originFile="NGS.A3SS.lst", targetFile="PB.A3SS.lst", outFile="Common.A3SS.lst")
#
#         makeLink("Common.SE.lst", "confident.SE.lst")
#         makeLink("Common.A5SS.lst", "confident.A5SS.lst")
#         makeLink("Common.A3SS.lst", "confident.A3SS.lst")
#
#         cmd = '''
#                 venn.R {NGS,PB}.SE.lst   -rD=90 -cP=180,0 -p=SE.pdf 2>/dev/null
#                 venn.R {NGS,PB}.A5SS.lst -rD=90 -cP=180,0 -p=A5SS.pdf 2>/dev/null
#                 venn.R {NGS,PB}.A3SS.lst -rD=90 -cP=180,0 -p=A3SS.pdf 2>/dev/null
#         '''
#         subprocess.call(cmd, shell=True)
#
#         filter(originFile="NGS.SE.lst", targetFile="PB.SE.lst", outFile="PB.specific.SE.lst")
#         filter(originFile="NGS.A5SS.lst", targetFile="PB.A5SS.lst", outFile="PB.specific.A5SS.lst")
#         filter(originFile="NGS.A3SS.lst", targetFile="PB.A3SS.lst", outFile="PB.specific.A3SS.lst")
#
#         cmd = "seDiagnose.pl -j {} > seDiagnose.tsv".format(pbSample.ngsJunctions)
#         subprocess.call(cmd, shell=True)
#         cmd = "anssDiagnose.pl -a 5 -j {} >a5ssDiagnose.tsv".format(pbSample.ngsJunctions)
#         subprocess.call(cmd, shell=True)
#         cmd = "anssDiagnose.pl -a 3 -j {} >a3ssDiagnose.tsv".format(pbSample.ngsJunctions)
#         subprocess.call(cmd, shell=True)
#     else:
#         makeLink("PB.SE.lst", "confident.SE.lst")
#         makeLink("PB.A5SS.lst", "confident.A5SS.lst")
#         makeLink("PB.A3SS.lst", "confident.A3SS.lst")
#
#     # cmd = '''
#     #         venn.R PB.IR.lst -p=IR.pdf 2>/dev/null
#     #         venn.R PB.APA.lst -p=APA.pdf 2>/dev/null
#     #         ln -sf PB.SE.lst confident.SE.lst
#     #         ln -sf PB.A5SS.lst confident.A5SS.lst
#     #         ln -sf PB.A3SS.lst confident.A3SS.lst
#     #     '''
#     # subprocess.call(cmd, shell=True)
#
#     # Classify APE into annotated or novel according to reference gene model
#     cmd = '''awk '{{print $0"\t"$12}}' {} | gpe2bed.pl -p >reference.bed12+'''.format(refParams.ref_gpe)
#     subprocess.call(cmd, shell=True, executable="/bin/bash")
#
#     findASmain(asType="IR", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="IR.reference.bed6+")
#     findASmain(asType="SE", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="SE.reference.bed6+")
#     findASmain(asType="A3SS", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="A3SS.reference.bed6+")
#     findASmain(asType="A5SS", gpeFile=refParams.ref_gpe, sourceFile="reference.bed12+", outFile="A5SS.reference.bed6+")
#     getAnnoASList("IR.reference.bed6+", "IR.reference.lst")
#     getAnnoASList("SE.reference.bed6+", "SE.reference.lst")
#     getAnnoASList("A3SS.reference.bed6+", "A3SS.reference.lst")
#     getAnnoASList("A5SS.reference.bed6+", "A5SS.reference.lst")
#     cmd = "paGroup.pl reference.bed12+ >paGroup.reference.tsv 2>paGroup.reference.bed6"
#     subprocess.call(cmd, shell=True)
#     cmd = "paCmp.pl -r paGroup.reference.tsv -a annoPA.tsv -n novelPA.tsv ../paGrouped.tsv >paGroup.novel.tsv 2>paGroup.anno.tsv"
#     subprocess.call(cmd, shell=True)
#
#     filter(originFile="IR.reference.lst", targetFile="PB.IR.lst", outFile="IR.anno.lst", mode="i")
#     filter(originFile="IR.reference.lst", targetFile="PB.IR.lst", outFile="IR.novel.lst", mode="e")
#     filter(originFile="SE.reference.lst", targetFile="confident.SE.lst", outFile="SE.anno.lst", mode="i")
#     filter(originFile="SE.reference.lst", targetFile="confident.SE.lst", outFile="SE.novel.lst", mode="e")
#     filter(originFile="A3SS.reference.lst", targetFile="confident.A3SS.lst", outFile="A3SS.anno.lst", mode="i")
#     filter(originFile="A3SS.reference.lst", targetFile="confident.A3SS.lst", outFile="A3SS.novel.lst", mode="e")
#     filter(originFile="A5SS.reference.lst", targetFile="confident.A5SS.lst", outFile="A5SS.anno.lst", mode="i")
#     filter(originFile="A5SS.reference.lst", targetFile="confident.A5SS.lst", outFile="A5SS.novel.lst", mode="e")
#
#     # Get Statistic of APE
#     getASstatistics(asType="IR", asFile="../PB/IR.bed6+", annoFile="IR.anno.lst", novelFile="IR.novel.lst",
#                     outFile="statistic.IR.tsv")
#     getASstatistics(asType="PA", annoFile="paGroup.anno.tsv", novelFile="paGroup.novel.tsv",
#                     outFile="statistic.APA.tsv")
#     getASstatistics(asType="SE", asFile="../PB/SE.confident.bed12+", annoFile="SE.anno.lst", novelFile="SE.novel.lst",
#                     outFile="statistic.SE.tsv")
#     getASstatistics(asType="A3SS", asFile="../PB/A3SS.confident.bed6+", annoFile="A3SS.anno.lst",
#                     novelFile="A3SS.novel.lst", outFile="statistic.A3SS.tsv")
#     getASstatistics(asType="A5SS", asFile="../PB/A5SS.confident.bed6+", annoFile="A5SS.anno.lst",
#                     novelFile="A5SS.novel.lst", outFile="statistic.A5SS.tsv")
#
#     # Event Decompose and Motif Evaluation
#     getSpliceSite(asType="IR", asFile="PB.IR.lst", outFile="IR.splicesite")
#     getSpliceSite(asType="IR", asFile="IR.anno.lst", outFile="IR.anno.splicesite")
#     getSpliceSite(asType="IR", asFile="IR.novel.lst", outFile="IR.novel.splicesite")
#     getSpliceSite(asType="SE", asFile="confident.SE.lst")
#     getSpliceSite(asType="A3SS", asFile="confident.A3SS.lst")
#     getSpliceSite(asType="A5SS", asFile="confident.A5SS.lst")
#     splicesite2figure(refParams=refParams)
#
#     getDist2TTS(refParams=refParams, paGroup="../PA/reads.paGrouped.bed6")
#
#     if pbSample.ngsProject != None:
#         cmd = '''
#             skyjoin <(cut -f4,5 ../PB/SE.bed12+) <(cut -f4,5 ../NGS/SE.bed12+) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=1-$3/1000;print $2,$3}' \
#             | correlation.R -x=PB -y=NGS -p=SE.cor.pdf >SE.cor 2>/dev/null
#             skyjoin <(cut -f4,5 ../PB/SE.bed12+) <(cut -f4,5 ../NGS/SE.bed12+) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=1-$3/1000;print "PB",$2,$1;print "NGS",$3,$1}' \
#             | line.R -b=V3 -s=0.01 -a=0.05 -x=Platform -y=PSI -m="R:$(head -n1 SE.cor | cut -f2), P:$(head -n1 SE.cor | cut -f3)" -p=SE.cor2.pdf 2>/dev/null
#             skyjoin <(cut -f4,5 ../PB/A5SS.bed6+) <(cut -f4,5 ../NGS/A5SS.{known,PB}.bed6+ | sort -u) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print $2,$3}' \
#             | correlation.R -x=PB -y=NGS -p=A5SS.cor.pdf >A5SS.cor 2>/dev/null
#             skyjoin <(cut -f4,5 ../PB/A5SS.bed6+) <(cut -f4,5 ../NGS/A5SS.{known,PB}.bed6+ | sort -u) 2>/dev/null \
#             | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print "PB",$2,$1;print "NGS",$3,$1}' \
#             | line.R -b=V3 -s=0.01 -a=0.05 -x=Platform -y=PSI -m="R:$(head -n1 A5SS.cor | cut -f2), P:$(head -n1 A5SS.cor | cut -f3)" -p=A5SS.cor2.pdf 2>/dev/null
#             skyjoin <(cut -f4,5 ../PB/A3SS.bed6+) <(cut -f4,5 ../NGS/A3SS.{known,PB}.bed6+ | sort -u) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print $2,$3}' \
#             | correlation.R -x=PB -y=NGS -p=A3SS.cor.pdf >A3SS.cor 2>/dev/null
#             skyjoin <(cut -f4,5 ../PB/A3SS.bed6+) <(cut -f4,5 ../NGS/A3SS.{known,PB}.bed6+ | sort -u) 2>/dev/null \
#             | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print "PB",$2,$1;print "NGS",$3,$1}' \
#             | line.R -b=V3 -s=0.01 -a=0.05 -x=Platform -y=PSI -m="R:$(head -n1 A3SS.cor | cut -f2), P:$(head -n1 A3SS.cor | cut -f3)" -p=A3SS.cor2.pdf 2>/dev/null
#         '''
#         subprocess.call(cmd, shell=True, executable="/bin/bash")
#     os.chdir(prevDir)
#     print str(datetime.datetime.now()) + " ASE Characterization for group {} done!".format(pbSample.group)

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    tgsAnalysisList = []
    project2as = {}
    pybedtools.set_tempdir(refParams.tmp_dir)
    try:
        pool = Pool(processes=len(projects))
        multiResults = []
        for project in projects:
            if project.tgsDataDir not in tgsAnalysisList:
                tgsAnalysisList.append(project.tgsDataDir)
                workDir = os.path.join(refParams.out_dir, project.projectName, project.group)
                if project.ngsProject:
                    project.ngsJunctions = os.path.join(workDir, "RNA-seq", "reassembly", "junctions.bed12")
                singleRunRes = pool.apply_async(identifyASE, (refParams, project, project2as))
                multiResults.append(singleRunRes)
                # characterizeASE(refParams=refParams, pbSample=project, project2as=project2as)
        for j in multiResults:
            j.wait()
    except Exception as e:
        print e

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)
