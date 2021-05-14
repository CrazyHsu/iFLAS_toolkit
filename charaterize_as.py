# from commonFuncs import *
from find_charaterize_as_functions import *
import pybedtools

def charaterize_as(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " ASE Characterization for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    characAsDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "as_events", "characterization")
    resolveDir(characAsDir)
    getAnnoASList("../ordinary_as/PB/SE.bed12+", "PB.SE.lst")
    getAnnoASList("../ordinary_as/PB/A5SS.bed6+", "PB.A5SS.lst")
    getAnnoASList("../ordinary_as/PB/A3SS.bed6+", "PB.A3SS.lst")
    getAnnoASList("../ordinary_as/PB/IR.bed6+", "PB.IR.lst")
    getAnnoASList("../pa/reads.paGrouped.tsv", "PB.APA.lst", PA=True)

    if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
        getAnnoASList("../ordinary_as/NGS/SE.bed12+", "NGS.SE.lst")
        getAnnoASList("../ordinary_as/NGS/A5SS.known.bed6+", "NGS.A5SS.lst")
        getAnnoASList("../ordinary_as/NGS/A5SS.PB.bed6+", "NGS.A5SS.lst", append=True, uniq=True)
        getAnnoASList("../ordinary_as/NGS/A3SS.known.bed6+", "NGS.A3SS.lst")
        getAnnoASList("../ordinary_as/NGS/A3SS.PB.bed6+", "NGS.A3SS.lst", append=True, uniq=True)

        filterFile(originFile="NGS.SE.lst", targetFile="PB.SE.lst", outFile="Common.SE.lst")
        filterFile(originFile="NGS.A5SS.lst", targetFile="PB.A5SS.lst", outFile="Common.A5SS.lst")
        filterFile(originFile="NGS.A3SS.lst", targetFile="PB.A3SS.lst", outFile="Common.A3SS.lst")

        makeLink("Common.SE.lst", "confident.SE.lst")
        makeLink("Common.A5SS.lst", "confident.A5SS.lst")
        makeLink("Common.A3SS.lst", "confident.A3SS.lst")

        cmd = '''
                venn.R {NGS,PB}.SE.lst   -rD=90 -cP=180,0 -p=SE.pdf 1>/dev/null 2>&1
                venn.R {NGS,PB}.A5SS.lst -rD=90 -cP=180,0 -p=A5SS.pdf 1>/dev/null 2>&1
                venn.R {NGS,PB}.A3SS.lst -rD=90 -cP=180,0 -p=A3SS.pdf 1>/dev/null 2>&1
        '''
        subprocess.call(cmd, shell=True)

        filterFile(originFile="NGS.SE.lst", targetFile="PB.SE.lst", outFile="PB.specific.SE.lst")
        filterFile(originFile="NGS.A5SS.lst", targetFile="PB.A5SS.lst", outFile="PB.specific.A5SS.lst")
        filterFile(originFile="NGS.A3SS.lst", targetFile="PB.A3SS.lst", outFile="PB.specific.A3SS.lst")

        cmd = "seDiagnose.pl -j {} PB.specific.SE.lst > seDiagnose.tsv".format(dataObj.ngs_junctions)
        subprocess.call(cmd, shell=True)
        cmd = "anssDiagnose.pl -a 5 -j {} PB.specific.A5SS.lst >a5ssDiagnose.tsv".format(dataObj.ngs_junctions)
        subprocess.call(cmd, shell=True)
        cmd = "anssDiagnose.pl -a 3 -j {} PB.specific.A3SS.lst >a3ssDiagnose.tsv".format(dataObj.ngs_junctions)
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
    findASmain(asType="IR", annoDict=annoDict, novelDict=novelDict, outFile="IR.reference.bed6+")
    findASmain(asType="SE", annoDict=annoDict, novelDict=novelDict, outFile="SE.reference.bed6+")
    findASmain(asType="A5SS", annoDict=annoDict, novelDict=novelDict, outFile="A5SS.reference.bed6+")
    findASmain(asType="A3SS", annoDict=annoDict, novelDict=novelDict, outFile="A3SS.reference.bed6+")

    getAnnoASList("IR.reference.bed6+", "IR.reference.lst")
    getAnnoASList("SE.reference.bed6+", "SE.reference.lst")
    getAnnoASList("A3SS.reference.bed6+", "A3SS.reference.lst")
    getAnnoASList("A5SS.reference.bed6+", "A5SS.reference.lst")
    cmd = "paGroup.pl reference.bed12+ >paGroup.reference.tsv 2>paGroup.reference.bed6"
    subprocess.call(cmd, shell=True)
    cmd = "paCmp.pl -r paGroup.reference.tsv -a annoPA.tsv -n novelPA.tsv ../pa/reads.paGrouped.tsv >paGroup.novel.tsv 2>paGroup.anno.tsv"
    subprocess.call(cmd, shell=True)

    filterFile(originFile="IR.reference.lst", targetFile="PB.IR.lst", outFile="IR.anno.lst", mode="i")
    filterFile(originFile="IR.reference.lst", targetFile="PB.IR.lst", outFile="IR.novel.lst", mode="e")
    filterFile(originFile="SE.reference.lst", targetFile="confident.SE.lst", outFile="SE.anno.lst", mode="i")
    filterFile(originFile="SE.reference.lst", targetFile="confident.SE.lst", outFile="SE.novel.lst", mode="e")
    filterFile(originFile="A3SS.reference.lst", targetFile="confident.A3SS.lst", outFile="A3SS.anno.lst", mode="i")
    filterFile(originFile="A3SS.reference.lst", targetFile="confident.A3SS.lst", outFile="A3SS.novel.lst", mode="e")
    filterFile(originFile="A5SS.reference.lst", targetFile="confident.A5SS.lst", outFile="A5SS.anno.lst", mode="i")
    filterFile(originFile="A5SS.reference.lst", targetFile="confident.A5SS.lst", outFile="A5SS.novel.lst", mode="e")

    # Get Statistic of APE
    getASstatistics(asType="IR", asFile="../ordinary_as/PB/IR.bed6+", annoFile="IR.anno.lst", novelFile="IR.novel.lst",
                    outFile="statistic.IR.tsv")
    getASstatistics(asType="PA", annoFile="paGroup.anno.tsv", novelFile="paGroup.novel.tsv",
                    outFile="statistic.APA.tsv")
    getASstatistics(asType="SE", asFile="../ordinary_as/PB/SE.confident.bed12+", annoFile="SE.anno.lst",
                    novelFile="SE.novel.lst", outFile="statistic.SE.tsv")
    getASstatistics(asType="A3SS", asFile="../ordinary_as/PB/A3SS.confident.bed6+", annoFile="A3SS.anno.lst",
                    novelFile="A3SS.novel.lst", outFile="statistic.A3SS.tsv")
    getASstatistics(asType="A5SS", asFile="../ordinary_as/PB/A5SS.confident.bed6+", annoFile="A5SS.anno.lst",
                    novelFile="A5SS.novel.lst", outFile="statistic.A5SS.tsv")

    # Event Decompose and Motif Evaluation
    getSpliceSite(asType="IR", asFile="PB.IR.lst", outFile="IR.splicesite")
    getSpliceSite(asType="IR", asFile="IR.anno.lst", outFile="IR.anno.splicesite")
    getSpliceSite(asType="IR", asFile="IR.novel.lst", outFile="IR.novel.splicesite")
    getSpliceSite(asType="SE", asFile="confident.SE.lst")
    getSpliceSite(asType="A3SS", asFile="confident.A3SS.lst")
    getSpliceSite(asType="A5SS", asFile="confident.A5SS.lst")
    splicesite2figure(refParams=refParams)

    getDist2TTS(refParams=refParams, paGroup="../pa/reads.paGrouped.bed6")

    # if dataObj.ngs_left_reads or dataObj.ngs_right_reads:
    #     cmd = '''
    #         skyjoin <(cut -f4,5 ../PB/SE.bed12+) <(cut -f4,5 ../NGS/SE.bed12+) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=1-$3/1000;print $2,$3}' \
    #         | correlation.R -x=PB -y=NGS -p=SE.cor.pdf >SE.cor 2>/dev/null
    #         skyjoin <(cut -f4,5 ../PB/SE.bed12+) <(cut -f4,5 ../NGS/SE.bed12+) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=1-$3/1000;print "PB",$2,$1;print "NGS",$3,$1}' \
    #         | line.R -b=V3 -s=0.01 -a=0.05 -x=Platform -y=PSI -m="R:$(head -n1 SE.cor | cut -f2), P:$(head -n1 SE.cor | cut -f3)" -p=SE.cor2.pdf 2>/dev/null
    #         skyjoin <(cut -f4,5 ../PB/A5SS.bed6+) <(cut -f4,5 ../NGS/A5SS.{known,PB}.bed6+ | sort -u) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print $2,$3}' \
    #         | correlation.R -x=PB -y=NGS -p=A5SS.cor.pdf >A5SS.cor 2>/dev/null
    #         skyjoin <(cut -f4,5 ../PB/A5SS.bed6+) <(cut -f4,5 ../NGS/A5SS.{known,PB}.bed6+ | sort -u) 2>/dev/null \
    #         | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print "PB",$2,$1;print "NGS",$3,$1}' \
    #         | line.R -b=V3 -s=0.01 -a=0.05 -x=Platform -y=PSI -m="R:$(head -n1 A5SS.cor | cut -f2), P:$(head -n1 A5SS.cor | cut -f3)" -p=A5SS.cor2.pdf 2>/dev/null
    #         skyjoin <(cut -f4,5 ../PB/A3SS.bed6+) <(cut -f4,5 ../NGS/A3SS.{known,PB}.bed6+ | sort -u) 2>/dev/null | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print $2,$3}' \
    #         | correlation.R -x=PB -y=NGS -p=A3SS.cor.pdf >A3SS.cor 2>/dev/null
    #         skyjoin <(cut -f4,5 ../PB/A3SS.bed6+) <(cut -f4,5 ../NGS/A3SS.{known,PB}.bed6+ | sort -u) 2>/dev/null \
    #         | awk 'BEGIN{OFS="\t"}{$2=$2/1000;$3=$3/1000;print "PB",$2,$1;print "NGS",$3,$1}' \
    #         | line.R -b=V3 -s=0.01 -a=0.05 -x=Platform -y=PSI -m="R:$(head -n1 A3SS.cor | cut -f2), P:$(head -n1 A3SS.cor | cut -f3)" -p=A3SS.cor2.pdf 2>/dev/null
    #     '''
    #     subprocess.call(cmd, shell=True, executable="/bin/bash")
    os.chdir(prevDir)
    print getCurrentTime() + " ASE Characterization for project {} entry {} done!".format(projectName, sampleName)
