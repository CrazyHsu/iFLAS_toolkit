import argparse, glob, os, pysam, itertools
import ConfigParser, StringIO, subprocess, datetime
from scipy.stats import chi2_contingency
import pandas as pd
from multiprocessing import Pool
# from Config import *
# from commonObjs import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

min_bq = 13
ploidy = 2

##############################################
def resolveDir(dirName):
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    os.chdir(dirName)

def getDictFromFile(myFile, sep="\t", inlineSep=None, keyCol=1, valueCol=None):
    with open(myFile) as f:
        myDict = {}
        for line in f.readlines():
            infoList = line.strip("\n").split(sep)
            key = infoList[keyCol-1]
            if valueCol:
                if inlineSep:
                    value = infoList[valueCol-1].split(inlineSep)
                else:
                    value = infoList[valueCol-1]
            else:
                value = infoList
            myDict[key] = value
        return myDict

def relationshipBetweenAlleleSpecifiAndAS(partial=False, nopartial=False, isoPairs=None):
    resultDict = {"partial": {}, "nopartial": {}}
    for isoIndex in range(len(isoPairs)):
        # partial
        if partial:
            partialHaplotype = "phased.partial.cleaned.human_readable.txt"
            partialDF = pd.read_csv(partialHaplotype, skiprows=[0], index_col=0, sep="\t")
            indexPairs = list(itertools.combinations(partialDF.index, 2))
            for pairIndex in range(len(indexPairs)):
                subPartialDF = partialDF.loc[indexPairs[pairIndex], isoPairs[isoIndex]]
                partialChi2Results = chi2_contingency(subPartialDF)
                partialChi2Pvalue = partialChi2Results[1]
                if partialChi2Pvalue <= 0.05:
                    combination = "_".join(map(str, ["partial", isoIndex, pairIndex]))
                    rows = subPartialDF.index
                    columns = subPartialDF.columns[pd.DataFrame(subPartialDF).idxmax()]
                    resultDict["partial"].update({combination: dict(zip(rows, columns))})

        # no-partial
        if nopartial:
            nopartialHaplotype = "phased.nopartial.cleaned.human_readable.txt"
            nopartialDF = pd.read_csv(nopartialHaplotype, skiprows=[0], index_col=0, sep="\t")
            indexPairs = list(itertools.combinations(nopartialDF.index, 2))
            for pairIndex in range(len(indexPairs)):
                subNopartialDF = nopartialDF.loc[indexPairs[pairIndex], isoPairs[isoIndex]]
                nopartialChi2Results = chi2_contingency(subNopartialDF)
                nopartialChi2Pvalue = nopartialChi2Results[1]
                if nopartialChi2Pvalue <= 0.05:
                    combination = "_".join(map(str, ["nopartial", isoIndex, pairIndex]))
                    rows = subNopartialDF.index
                    columns = subNopartialDF.columns[pd.DataFrame(subNopartialDF).idxmax()]
                    resultDict["nopartial"].update({combination: dict(zip(rows, columns))})

        if not partial and nopartial:
            raise Exception("You must specify either partial or nopartial file")
    return resultDict

def getASpairedIsoforms(asFile, collapsedGroupFile, asType="SE", filterByCount=0):
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=1)
    asPairs = {}
    with open(asFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            if asType == "SE":
                inclusionIsos = records[15].split(",")
                exclusionIsos = records[17].split(",")
            elif asType == "PA":
                inclusionIsos = []
                exclusionIsos = []
            else:
                inclusionIsos = records[7].split(",")
                exclusionIsos = records[9].split(",")

            for item in itertools.product(inclusionIsos, exclusionIsos):
                inclusionReads = collapsedTrans2reads[item[0]]
                exclusionReads = collapsedTrans2reads[item[1]]
                if filterByCount:
                    if len(inclusionReads) < filterByCount or len(exclusionReads) < filterByCount:
                        continue
                '''the isoform PB.x.x split by "." to determine gene PB.x'''
                inclusionGene = ".".join(item[0].split(".")[:2])
                exclusionGene = ".".join(item[1].split(".")[:2])
                isoPair = [item[0], item[1]]
                if inclusionGene == exclusionGene:
                    if inclusionGene not in asPairs:
                        asPairs[inclusionGene] = [isoPair]
                    else:
                        asPairs[inclusionGene].append(isoPair)
        return asPairs

def makeAbundanceFile(groupFile, outFile=None):
    with open(groupFile) as f:
        cidInfo = {}
        for line in f.readlines():
            pbid, members = line.strip().split('\t')
            for cid in members.split(','):
                cidInfo[cid] = pbid
        if outFile:
            out = open(outFile, "w")
            print >> out, "id\tlength\tis_fl\tstat\tpbid"
            for i in cidInfo:
                print >> out, "\t".join([i, "NA", "Y", "unique", cidInfo[i]])
            out.close()
        else:
            print "id\tlength\tis_fl\tstat\tpbid"
            for i in cidInfo:
                print "\t".join([i, "NA", "Y", "unique", cidInfo[i]])

def runPhaser(targetDir, asPairs):
    curDir = os.getcwd()
    os.chdir(targetDir)
    configFile = "config"
    with open(configFile) as f:
        config_string = StringIO.StringIO("[dummy_section]\n" + f.read())
    config_parser = ConfigParser.RawConfigParser()
    config_parser.readfp(config_string)
    strand = config_parser.get("dummy_section", "ref_strand")

    cmd = "minimap2 -ax splice fake.fasta ccs.fastq >ccs.sam 2>/dev/null"
    subprocess.call(cmd, shell=True)

    pysam.sort("-o", "ccs.sorted.bam", "-@", str(5), "ccs.sam", catch_stdout=False)
    pileupRes = pysam.mpileup("--min-BQ", str(min_bq), "-f", "fake.fasta", "-s", "ccs.sorted.bam")
    pileupResOut = open("ccs.mpileup", "w")
    pileupResOut.write(pileupRes)
    pileupResOut.close()

    cmd = "run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --strand {} -o phased.nopartial -n {} 1>nopartial.stdout.txt".format(
        strand, ploidy)
    subprocess.call(cmd, shell=True)
    cmd = "run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --partial_ok --strand {} -o phased.partial -n {} 1>partial.stdout.txt".format(
        strand, ploidy)
    subprocess.call(cmd, shell=True)

    os.chdir(curDir)


def alleleSpecific_test():
    collapsedGroupFile = "../filtration/collapse/tofu.collapsed.group.txt"
    readStatFile = "tofu.collapsed.read_stat.txt"
    makeAbundanceFile(collapsedGroupFile, outFile=readStatFile)

    seFile = "../ASE/PB/SE.confident.bed12+"
    irFile = "../ASE/PB/IR.confident.bed6+"
    a5ssFile = "../ASE/PB/A5SS.confident.bed6+"
    a3ssFile = "../ASE/PB/A3SS.confident.bed6+"

    seAsPairs = getASpairedIsoforms(seFile, collapsedGroupFile, asType="SE")
    irAsPairs = getASpairedIsoforms(irFile, collapsedGroupFile, asType="IR")
    a5ssAsPairs = getASpairedIsoforms(a5ssFile, collapsedGroupFile, asType="A5SS")
    a3ssAsPairs = getASpairedIsoforms(a3ssFile, collapsedGroupFile, asType="A3SS")

    asPairs = {"SE": seAsPairs, "IR": irAsPairs, "A5SS": a5ssAsPairs, "A3SS": a3ssAsPairs}

    ref_genome = "/data/Zea_mays.AGPv3.25/Zea_mays.AGPv3.25.dna.genome.fa"
    flncFq = "../isoseq3/flnc.fq"
    collapsedGff = "../filtration/collapse/tofu.collapsed.gff"
    cmd = "select_loci_to_phase.py {} {} {} {} -c 40 1>/dev/null 2>&1".format(ref_genome, flncFq,
                                                                              collapsedGff, readStatFile)
    subprocess.call(cmd, shell=True)

    lociDir = glob.glob("by_loci/*size*")
    runPhaser(lociDir[0], asPairs)

def alleleSpecific(refParams=None, paramObj=None):
    print str(datetime.datetime.now()) + " Analysis the relationship between allele specific and AS events for project {} group {}...".format(
        paramObj.projectName, paramObj.group)
    # prevDir = os.getcwd()
    groupDir = os.path.join(refParams.out_dir, paramObj.projectName, paramObj.group)
    resolveDir(os.path.join(groupDir, "alleleSpecificRelatedAS"))

    collapsedGroupFile = os.path.join(groupDir, "filtration", "collapse", "tofu.collapsed.group.txt")
    readStatFile = "tofu.collapsed.read_stat.txt"
    makeAbundanceFile(collapsedGroupFile, outFile=readStatFile)

    # select_loci_to_phase = os.path.join(cDNA_Cupcake_dir, "phasing", "utils", "select_loci_to_phase.py")
    flncFq = os.path.join(groupDir, "isoseq3", "flnc.fq")
    collapsedGff = os.path.join(groupDir, "filtration", "collapse", "tofu.collapsed.gff")

    seFile = os.path.join(groupDir, "ASE", "PB", "SE.confident.bed12+")
    irFile = os.path.join(groupDir, "ASE", "PB", "IR.confident.bed6+")
    a5ssFile = os.path.join(groupDir, "ASE", "PB", "A5SS.confident.bed6+")
    a3ssFile = os.path.join(groupDir, "ASE", "PB", "A3SS.confident.bed6+")

    seAsPairs = getASpairedIsoforms(seFile, collapsedGroupFile, asType="SE")
    irAsPairs = getASpairedIsoforms(irFile, collapsedGroupFile, asType="IR")
    a5ssAsPairs = getASpairedIsoforms(a5ssFile, collapsedGroupFile, asType="A5SS")
    a3ssAsPairs = getASpairedIsoforms(a3ssFile, collapsedGroupFile, asType="A3SS")

    asPairs = {"SE": seAsPairs, "IR": irAsPairs, "A5SS": a5ssAsPairs, "A3SS": a3ssAsPairs}

    # generate "by_loci" directory
    # removeFiles()
    cmd = "select_loci_to_phase.py {} {} {} {} -c 40 1>/dev/null 2>&1".format(refParams.ref_genome, flncFq, collapsedGff, readStatFile)
    # subprocess.call(cmd, shell=True)

    lociDir = glob.glob("by_loci/*size*")
    # multiResults = []
    # pool = Pool(processes=refParams.threads)
    # for i in lociDir:
    #     singleRunRes = pool.apply_async(runPhaser, (i, refParams, asPairs))
    #     multiResults.append(singleRunRes)
    # for j in multiResults:
    #     j.wait()

    resultDict = {}
    for i in lociDir:
        os.chdir(i)
        fakeGene = os.path.basename(i).split("_")[0]
        resultDict[fakeGene] = {}
        for asType in asPairs:
            if fakeGene not in asPairs[asType]:
                continue
            isoPairs = asPairs[asType][fakeGene]
            sigRelatedDict = relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=True, isoPairs=isoPairs)
            resultDict[fakeGene].update({asType: sigRelatedDict})

    '''the structure of resultDict is resultDict = {fakeGene: {asType: {partialCategory: {combinationName: {haplotype1: PB.X.1, haplotype2: PB.X.2, ...}}}}}'''
    partialOut = open("partialAsRelatedHaplotype.txt", "w")
    nopartialOut = open("nopartialAsRelatedHaplotype.txt", "w")
    for g in resultDict:
        for a in resultDict[g]:
            for p in resultDict[g][a]:
                for c in resultDict[g][a][p]:
                    for h in resultDict[g][a][p][c]:
                        if p == "partial":
                            print >> partialOut, "\t".join([g, a, h, resultDict[g][a][p][c][h], c])
                        else:
                            print >> nopartialOut, "\t".join([g, a, h, resultDict[g][a][p][c][h], c])
    partialOut.close()
    nopartialOut.close()
    # relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=True, asPairs=seAsPairs)
    # relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=True, asPairs=irAsPairs)
    # relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=True, asPairs=a5ssAsPairs)
    # relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=True, asPairs=a3ssAsPairs)


    print str(datetime.datetime.now()) + " Analysis the relationship between allele specific and AS events for project {} group {} done!".format(
        paramObj.projectName, paramObj.group)

# def main(defaultCfg=None):
#     researchCfg = defaultCfg.researchCfg
#     refParams = defaultCfg.refParams
#     projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
#     projects = projectCfg.projects
#     projectsDict = projectCfg.projectsDict
#     initSysSetting(refParams=refParams, projectsDict=projectsDict)
#     tgsAnalysisList = []
#     poolNum = sum([len(projectsDict[i].keys()) for i in projectsDict])
#     try:
#         pool = MyPool(processes=poolNum)
#         multiResults = []
#         for project in projects:
#             if project.tgsDataDir not in tgsAnalysisList:
#                 tgsAnalysisList.append(project.tgsDataDir)
#                 singleRunRes = pool.apply_async(alleleSpecific, (refParams, project))
#                 multiResults.append(singleRunRes)
#         for j in multiResults:
#             j.wait()
#     except Exception as e:
#         print e

if __name__ == '__main__':
    alleleSpecific_test()
    # defaultCfg = Config1(args.default_cfg)
    # main(defaultCfg=defaultCfg)