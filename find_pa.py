# from commonFuncs import *
from find_charaterize_as_functions import *
import pandas as pd

def find_pa(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Polyadenylation analysis for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    paDir = os.path.join(baseDir, "as_events", "pa")
    resolveDir(paDir)
    makeLink(os.path.join(baseDir, "collapse", "reads.assigned.unambi.bed12+"), "reads.assigned.unambi.bed12+")
    makeLink(os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt"), "tofu.collapsed.group.txt")
    makeLink(os.path.join(baseDir, "collapse", "isoformGrouped.bed12+"), "isoformGrouped.bed12+")

    with open("reads.assigned.unambi.bed12+") as f:
        cleavageList = {}
        count = 0
        cleavageOut = open("cleavage.bed6", "w")
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            chrom, start, end, strand = lineInfo[0], int(lineInfo[1]), int(lineInfo[2]), lineInfo[5]
            if strand == "+":
                keyStr = "{}:{}-{}:{}".format(chrom, end-1, end, strand)
                if keyStr not in cleavageList:
                    count += 1
                    cleavageList.update({keyStr: 1})
                    print >>cleavageOut, "\t".join(map(str, [chrom, end - 1, end, "ClevageSite{}".format(count), 0, strand]))
            elif strand == "-":
                keyStr = "{}:{}-{}:{}".format(chrom, start, start + 1, strand)
                if keyStr not in cleavageList:
                    count += 1
                    cleavageList.update({keyStr: 1})
                    print >>cleavageOut, "\t".join(map(str, [chrom, start, start + 1, "ClevageSite{}".format(count), 0, strand]))
        cleavageOut.close()

    cmd = '''bedClosest.pl cleavage.bed6 | tee adjacentCleavage.tsv | cut -f13 | distrCurve.R -xl=2 -d -x='Binned Distance between Adjacent Cleavage Site' -p=adjacentCleavageDistr.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)

    # pa(polish_flnc_cluster="clustered_unclustered.merged_report.csv", bed12="target.transcript.correlated.flnc.sorted.bed12+",
    #    tofu_group="tofu.collapsed.group.txt", threads=refParams.threads, out="paCluster.bed8+")
    getPaCluster(readsBed="reads.assigned.unambi.bed12+", tofuGroup="tofu.collapsed.group.txt", threads=dataObj.single_run_threads, paClusterOut="paCluster.bed8+")
    cmd = '''awk '$5>1{print $3-$2}' paCluster.bed8+ |
             distrCurve.R -d -x='Cluster Size (limited in 1-100)' -y='Density' -m='Distribution of Cluster Size' -x1=0 -x2=100 -b=1 -p=paClusterSize.pdf 2>/dev/null
    '''
    subprocess.call(cmd, shell=True)
    paBed6 = open("PA.bed6+", "w")
    with open("paCluster.bed8+") as f:
        for line in f.readlines():
            lineInfo = line.strip("\n").split("\t")
            print >> paBed6, "\t".join(map(str, [lineInfo[0], lineInfo[6], lineInfo[7]] + lineInfo[3:6] + lineInfo[8:]))
    paBed6.close()
    cmd = "paGroup.pl isoformGrouped.bed12+ >isoform.paGrouped.tsv 2>isoform.paGrouped.bed6"
    subprocess.call(cmd, shell=True)
    cmd = '''3endRevise.pl -p PA.bed6+ <(cut -f 1-12,15 reads.assigned.unambi.bed12+) | tee reads.3endRevised.bed12+ | paGroup.pl >reads.paGrouped.tsv 2>reads.paGrouped.bed6'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    cmd = "PAClassbyRead.pl -a reads.assigned.unambi.bed12+ <(cut -f1-8 paCluster.bed8+) >paCluster.type.bed8+ 2>singleExonReadWithExonInMEread.bed12+"
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    deSingleExon = open("paCluster.deSingleExonRead.bed8+", "w")
    singleExon = open("paCluster.singleExonRead.bed8+", "w")
    with open("paCluster.type.bed8+") as f:
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            if lineInfo[8] != "SE":
                print >> deSingleExon, line.strip("\n")
            elif lineInfo[8] == "SE":
                print >> singleExon, line.strip("\n")
    deSingleExon.close()
    singleExon.close()

    resolveDir("motif")
    motifAroundPA("../PA.bed6+", up1=100, down1=100, up2=50, down2=0, refFasta=refParams.ref_genome, chrLenFile=refParams.ref_size)

    cmd = '''lines.R -w=15 -y=Frequency -x='Distance Relative to PA' -m='Distribution of Nucleotide Frequency around PA' -p=nucleotide.pdf *.nucleotide 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    cmd = '''lines.R -p=PAS.1.pdf {AATAAA,AAATAA,ATAAAA,ATTAAA,ATAAAT,TAATAA}.PAS -x1=-50 -x2=0 -w=15 2>/dev/null &&
             lines.R -p=PAS.2.pdf {ATAAAG,AAAATA,CAATAA,ATAAAC,AAAAAA,AAAAAG}.PAS -x1=-50 -x2=0 -w=15 2>/dev/null
    '''
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    otherMotifFileList = ["ATAAAC.PAS", "ATAAAG.PAS", "CAATAA.PAS", "AAAAAA.PAS", "AAAAAG.PAS", "AAAATA.PAS"]
    motifPercentDF = [pd.read_csv(i, sep="\t", index_col=0, header=None) for i in otherMotifFileList]
    motifPercentSummedDF = pd.concat(motifPercentDF, axis=1).sum(axis=1)
    motifPercentSummedDF.to_frame().to_csv("Other.PAS", sep="\t", header=None)

    cmd = '''lines.R -p=PAS.pdf {Other,AATAAA,AAATAA,ATAAAA,ATTAAA,ATAAAT,TAATAA}.PAS -x1=-50 -x2=0 -w=15 2>/dev/null'''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Polyadenylation analysis for project {} entry {} done!".format(projectName, sampleName)
