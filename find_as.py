#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: find_as.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:16:33
Last modified: 2021-04-29 16:16:33
'''
from commonFuncs import *
# from commonObjs import *
from find_as_functions import *
import pybedtools



def find_perl():
    cmd = '''
    getSingleIr.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/IR.bed6+
geneWithIR=$(cut -f7 PB/IR.bed6+ | sort -u | wc -l)
deSingleExonReadCount=$(cut -f13 deSingleExonRead.bed12+ | sort -u | wc -l)
geneWithIrRatio=$(awk 'BEGIN{printf "%.4f", '$geneWithIR/$deSingleExonReadCount'}')
getSE.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/SE.bed12+
getA5SS.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/A5SS.bed6+
getA3SS.pl -g ../../raw.gpe deSingleExonRead.bed12+ >PB/A3SS.bed6+
awk '$9>1 && $11>1 && $5>=100 && $5<=900{print $3-$2}' PB/A5SS.bed6+ | hist.R -p=PB/A5SS.size.pdf 2>/dev/null
awk '$3-$2>=20 && $9>1 && $11>1{print $5/1000}' PB/A5SS.bed6+ | hist.R -x='Incusion Ratio' -p=PB/A5SS.InclusionRatio.pdf 2>/dev/null
ln -sf SE.bed12+ PB/SE.confident.bed12+
ln -sf A5SS.bed6+ PB/A5SS.confident.bed6+
ln -sf A3SS.bed6+ PB/A3SS.confident.bed6+
awk '{print $3-$2}' PB/A5SS.confident.bed6+ | hist.R -p=PB/A5SS.size.pdf 2>/dev/null
awk '$9>1 && $11>1{print $5/1000}' PB/A5SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=PB/A5SS.InclusionRatio.pdf 2>/dev/null
awk '{print $3-$2}' PB/A3SS.confident.bed6+ | hist.R -p=PB/A3SS.size.pdf 2>/dev/null
awk '$9>1 && $11>1{print $5/1000}' PB/A3SS.confident.bed6+ | hist.R -x='Inclusion Ratio' -p=PB/A3SS.InclusionRatio.pdf 2>/dev/null
echo -e "\n[INFO] $(date) Alternative splicing events identifying done"

    '''

def find_python(isoformBed, tofuGroupFile, dataObj=None, refParams=None):
    if not os.path.exists("PB"):
        os.makedirs("PB")
    transBedList = GenePredObj(refParams.ref_gpe, bincolumn=False).toBed(gene=True)
    transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)
    readsBedList = []
    isoform2reads = {}
    with open(tofuGroupFile) as f:
        for line in f.readlines():
            isoform, reads = line.strip("\n").split("\t")
            isoform2reads[isoform] = reads.split(",")
    with open(isoformBed) as f:
        for line in f.readlines():
            readStruc = Bed12(line)
            if len(isoform2reads[readStruc.name]) < 2: continue
            readsBedList.append("\t".join(readStruc.record[:12]))

    readsBedObj = pybedtools.BedTool("\n".join(readsBedList), from_string=True)

    annoBedRes = readsBedObj.intersect(transBedObj, wa=True, wb=True, s=True)
    novelBedRes = readsBedObj.intersect(transBedObj, v=True, s=True)
    annoDict, novelDict = getAS(annoBedRes, novelBedRes, offset=0)
    findASmain(asType="IR", annoDict=annoDict, novelDict=novelDict, outFile="PB/IR.bed6+",
                   isoform2reads=isoform2reads)
    findASmain(asType="SE", annoDict=annoDict, novelDict=novelDict, outFile="PB/SE.bed12+",
                   isoform2reads=isoform2reads)
    findASmain(asType="A5SS", annoDict=annoDict, novelDict=novelDict, outFile="PB/A5SS.bed6+",
                   isoform2reads=isoform2reads)
    findASmain(asType="A3SS", annoDict=annoDict, novelDict=novelDict, outFile="PB/A3SS.bed6+",
                   isoform2reads=isoform2reads)

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

    if dataObj.ngsLeftReads or dataObj.ngsRightReads:
        if not os.path.exists("NGS"):
            os.makedirs("NGS")
        try:
            from scanESbyNGS import scanEsByNGS
            scanEsByNGS(fref=refParams.ref_gpe, has_bin=False, fjunc=dataObj.ngsJunctions,
                        fpb="isoformGrouped.bed12+", outFile="SEsupportByNGS.bed12+")
            assignGeneNameToSE(fref=refParams.ref_gpe, has_bin=False, fes="SEsupportByNGS.bed12+",
                               outFile="NGS/SE.bed12+", errorOutFile="NGS/novelSE2gene.log")
        except Exception as e:
            print e
        filterFile(originFile="NGS/SE.bed12+", targetFile="PB/SE.bed12+", originField=4, targetField=4,
               outFile="PB/SE.NGS.bed12+")
        filterFile(originFile="PB/SE.bed12+", targetFile="NGS/SE.bed12+", originField=4, targetField=4,
               outFile="NGS/SE.PB.bed12+")
        cmd = "juncAssign.pl -g {} {} | tee junction.assigned.for_IR.bed12+ | cut -f 1-12,21-22 >junction.assigned.bed12+".format(
            refParams.ref_gpe, dataObj.ngsJunctions)
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
    # os.chdir(prevDir)
    # characterizeASE(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)


def find_as(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Alternative splicing events identifying for project {} sample {}...".format(
        projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dataObj.out_dir, projectName, sampleName)
    logDir = os.path.join(baseDir, "log")

    workDir = os.path.join(baseDir, "as_events")
    resolveDir(workDir)
    isoformBed = os.path.join(baseDir, "collapse", "isoformGrouped.bed12+")
    tofuGroupFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    if dataObj.ngsJunctions == None:
        dataObj.ngsJunctions = os.path.join(baseDir, "rna-seq", "reassembly", "junctions.bed")

    if not os.path.exists("PB"):
        os.makedirs("PB")
    find_python(isoformBed, tofuGroupFile, dataObj=dataObj, refParams=refParams)
    # transBedList = GenePredObj(refParams.ref_gpe, bincolumn=False).toBed(gene=True)
    # transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)
    # readsBedList = []
    # isoform2reads = {}
    # with open("tofu.collapsed.group.txt") as f:
    #     for line in f.readlines():
    #         isoform, reads = line.strip("\n").split("\t")
    #         isoform2reads[isoform] = reads.split(",")
    # with open("isoformGrouped.bed12+") as f:
    #     for line in f.readlines():
    #         readStruc = Bed12(line)
    #         if len(isoform2reads[readStruc.name]) < 2: continue
    #         readsBedList.append("\t".join(readStruc.record[:12]))

    os.chdir(prevDir)
    print getCurrentTime() + " Alternative splicing events identifying for project {} sample {}...".format(
        projectName, sampleName)
