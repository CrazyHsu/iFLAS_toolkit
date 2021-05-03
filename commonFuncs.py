import os, re, datetime, psutil, subprocess, pysam
from Bio import SeqIO

def batchCreateDir(dirList):
    for i in dirList:
        if not os.path.exists(i):
            os.makedirs(i)

def checkFast5Files(fast5Dir):
    if os.path.isfile(fast5Dir):
        if fast5Dir.endswith("fq") or fast5Dir.endswith("fastq"):
            return "fastq"
    else:
        fileCount = len([x for x in os.listdir(fast5Dir) if not x.endswith(".txt")])
        fast5FileCount = 0
        fqFileCount = 0
        for i in os.listdir(fast5Dir):
            if i.endswith("fast5"):
                fast5FileCount += 1
            if i.endswith("fq") or i.endswith("fastq"):
                fqFileCount += 1
        if fast5FileCount == fileCount:
            return "fast5"
        if fqFileCount == fileCount:
            return "fastq"
        if fast5FileCount != fileCount or fqFileCount != fileCount:
            raise Exception("Please check the {} directory which is not a pure fast5 or fastq-containing directory".format(fast5Dir))

def checkHisat2IndexExist(dir):
    if len(os.listdir(dir)) == 0:
        return False
    else:
        ht2Flag = 0
        for i in os.listdir(dir):
            if i.endswith("exon") or i.endswith("ss") or i.endswith("log"):
                continue
            elif i.endswith("ht2"):
                ht2Flag = 1
                continue
            else:
                return False
        if ht2Flag:
            return True
        else:
            return False

def makeHisat2Index(refParams=None, hisat2indexDir=None, indexPrefix=None, optionTools=None):
    # hisat2indexDir = os.path.join(outDir, "hisat2indexDir")
    if not os.path.exists(hisat2indexDir):
        os.makedirs(hisat2indexDir)

    if not checkHisat2IndexExist(hisat2indexDir):
        print str(datetime.datetime.now()) + " Start Hisat2 indexing..."
        batchCreateDir([hisat2indexDir])
        # gtfPrefix = os.path.splitext(os.path.basename(refParams.ref_gtf))[0]
        refAnnoGTF = refParams.ref_gtf
        refAnnoSS = "{}/{}.ss".format(hisat2indexDir, indexPrefix)
        cmd = "hisat2_extract_splice_sites.py {} >{}".format(refAnnoGTF, refAnnoSS)
        subprocess.call(cmd, shell=True)
        refAnnoExon = "{}/{}.exon".format(hisat2indexDir, indexPrefix)
        cmd = "hisat2_extract_exons.py {} >{}".format(refAnnoGTF, refAnnoExon)
        subprocess.call(cmd, shell=True)
        cmd = "hisat2-build -p {} --ss {} --exon {} {} {}/{} 1>{}/hisat2build.log 2>&1".format(optionTools.threads,
                                                                                               refAnnoSS, refAnnoExon,
                                                                                               refParams.ref_genome,
                                                                                               hisat2indexDir,
                                                                                               indexPrefix, hisat2indexDir)
        subprocess.call(cmd, shell=True)
        print str(datetime.datetime.now()) + " End Hisat2 indexing!"

def getSubSamByName(samFile, nameList, isBam=False, nameListIsFile=False, outPrefix=None, sort=True, threads=4):
    n = nameList
    if nameListIsFile:
        with open(nameList, 'r') as infile:
            n = infile.read().splitlines()
    if isBam:
        samHandle = pysam.AlignmentFile(samFile, 'rb', check_sq=False)
        suffix = ".bam"
    else:
        tmpBam = "tmp.bam"
        pysam.view("-o", tmpBam, "--output-fmt", "BAM", samFile, catch_stdout=False)
        samHandle = pysam.AlignmentFile(tmpBam, 'rb', check_sq=False)
        suffix = ".sam"
    name_indexed = pysam.IndexedReads(samHandle)
    name_indexed.build()
    header = samHandle.header.copy()
    if isBam:
        out = pysam.Samfile(outPrefix + suffix, 'wb', header=header)
    else:
        out = pysam.Samfile(outPrefix + suffix, 'w', header=header)
    for name in n:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                out.write(x)
    out.close()
    if not isBam:
        os.remove("tmp.bam")
    if sort:
        if isBam:
            pysam.sort("-o", outPrefix + ".sorted" + suffix, "-@", str(threads), outPrefix + suffix, catch_stdout=False)
        else:
            pysam.sort("-o", outPrefix + ".sorted" + suffix, "-@", str(threads), "--output-fmt", "SAM", outPrefix + suffix, catch_stdout=False)
        os.remove(outPrefix + suffix)
        return outPrefix + ".sorted" + suffix
    else:
        return outPrefix + suffix


def getFxSequenceId(fxFile, isFa=False, isFq=False, outFile=None):
    recordId = []
    if isFa:
        for rec in SeqIO.parse(fxFile, "fasta"):
            recordId.append(rec.id)
        if outFile:
            with open(outFile, 'w') as f:
                for item in recordId:
                    print >> f, item
        else:
            return recordId
    if isFq:
        for rec in SeqIO.parse(fxFile, "fastq"):
            recordId.append(rec.id)
        if outFile:
            with open(outFile, 'w') as f:
                for item in recordId:
                    print >> f, item
        else:
            return recordId

def getCurrentTime():
    return str(datetime.datetime.now())

def initRefSetting(refParams=None, dirSpec=None):
    originDir = os.getcwd()
    if refParams.ref_genome == None:
        raise Exception("You must specify a genome file!")
    refParams.ref_genome = os.path.abspath(os.path.join(originDir, refParams.ref_genome))
    validateFile(refParams.ref_genome)

    if refParams.ref_gtf == None and refParams.ref_gff == None:
        raise Exception("You must specify a GTF or GFF3 file!")
    elif refParams.ref_gtf == None and refParams.ref_gff != None:
        refParams.ref_gff = os.path.abspath(os.path.join(originDir, refParams.ref_gff))
        validateFile(refParams.ref_gff)
        gffDir = os.path.dirname(refParams.ref_gff)
        ref_gtf = os.path.abspath(os.path.join(gffDir, "iFLAS.gtf"))
        cmd = "gffread -T -o {} {}".format(ref_gtf, refParams.ref_gff)
        subprocess.call(cmd, shell=True)
        refParams.ref_gtf = ref_gtf

    gtfDir = os.path.dirname(refParams.ref_gtf)
    if refParams.ref_size == None:
        refParams.ref_size = os.path.join(gtfDir, "iFLAS.len.txt")
        sizeOut = open(refParams.ref_size)
        for seqRecord in SeqIO.parse(refParams.ref_genome, "fasta"):
            print >> sizeOut, "\t".join([seqRecord.id, seqRecord.seq])
        sizeOut.close()
    else:
        refParams.ref_size = os.path.join(originDir, refParams.ref_size)

    if refParams.ref_gpe == None:
        refParams.ref_gpe = os.path.join(gtfDir, "iFLAS.gpe")
        cmd = "gtfToGenePred {} {} -genePredExt".format(refParams.ref_gtf, refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
    else:
        refParams.ref_gpe = os.path.abspath(os.path.join(originDir, refParams.ref_gpe))

    if refParams.ref_bed == None:
        refParams.ref_bed = os.path.join(gtfDir, "iFLAS.annotation.bed")
        cmd = "gpe2bed.pl -g {} > {}".format(refParams.ref_gpe, refParams.ref_bed)
        subprocess.call(cmd, shell=True)
    else:
        refParams.ref_bed = os.path.abspath(os.path.join(originDir, refParams.ref_bed))
    refParams.ref_mm2_index = os.path.abspath(os.path.join(originDir, refParams.ref_mm2_index))

    if dirSpec.tmp_dir == None:
        dirSpec.tmp_dir = "/tmp"
    else:
        if not os.path.exists(dirSpec.tmp_dir):
            batchCreateDir([dirSpec.tmp_dir])
    dirSpec.out_dir = os.path.abspath(os.path.join(originDir, dirSpec.out_dir))
    if not os.path.exists(dirSpec.out_dir):
        batchCreateDir([dirSpec.out_dir])

def initSysResourceSetting(optionTools=None):
    maxMem = psutil.virtual_memory().total >> 20
    if str(optionTools.memory).endswith("M"):
        if int(optionTools.memory[0:-1]) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    elif str(optionTools.memory).endswith("G"):
        optionTools.memory = str(int(optionTools.memory[0:-1]) * 1024) + "M"
        if int(optionTools.memory[0:-1]) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    elif optionTools.memory.isdigit():
        if int(optionTools.memory) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    else:
        raise Exception("Please check the format of memory you specified in optionTool section!")

    maxThreads = psutil.cpu_count()
    if str(optionTools.threads).isdigit():
        if int(optionTools.threads) > maxThreads * 0.75:
            optionTools.threads = int(maxThreads * 0.75)
        else:
            optionTools.threads = int(optionTools.threads)
    else:
        raise Exception("The number of threads should be digital, please check it!")


def resolveDir(dirName, chdir=True):
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    if chdir:
        os.chdir(dirName)

def removeFiles(myDir, fileList):
    for f in fileList:
        os.remove(os.path.join(myDir, f.strip("\n")))
        
def renameNGSdata2fastp(dataObj=None):
    if dataObj.ngsPaired == "paired":
        leftReadsRepeats = [i.strip() for i in dataObj.ngsLeftReads.split(";")]
        rightReadsRepeats = [i.strip() for i in dataObj.ngsRightReads.split(";")]
        if len(leftReadsRepeats) != len(rightReadsRepeats):
            raise Exception("The repeats of your NGS data not match between your left reads and right reads")
        else:
            newLeftReadsRepeats = []
            newRightReadsRepeats = []
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                rightReads = rightReadsRepeats[i].split(",")
                if len(leftReads) != len(rightReads):
                    raise Exception("please input the right paired NGS reads for the analysis")
                else:
                    newLeftReads = []
                    newRightReads = []
                    for j in range(len(leftReads)):
                        leftReadsDir = os.path.dirname(leftReads[j])
                        rightReadsDir = os.path.dirname(rightReads[j])
                        leftReadsBase = os.path.basename(leftReads[j])
                        rightReadsBase = os.path.basename(rightReads[j])
                        newLeft = os.path.join(leftReadsDir, "fastp.{}".format(leftReadsBase))
                        newRight = os.path.join(rightReadsDir, "fastp.{}".format(rightReadsBase))
                        newLeftReads.append(newLeft)
                        newRightReads.append(newRight)
                    newLeftReadsRepeats.append(",".join(newLeftReads))
                    newRightReadsRepeats.append(",".join(newRightReads))
            dataObj.ngsLeftReads = ";".join(newLeftReadsRepeats)
            dataObj.ngsRightReads = ";".join(newRightReadsRepeats)
    else:
        if dataObj.ngsLeftReads and dataObj.ngsRightReads == None:
            leftReadsRepeats = [i.strip() for i in dataObj.ngsLeftReads.split(";")]
            newLeftReadsRepeats = []
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                newLeftReads = []
                for j in range(len(leftReads)):
                    leftReadsDir = os.path.dirname(leftReads[j])
                    leftReadsBase = os.path.basename(leftReads[j])
                    newLeft = os.path.join(leftReadsDir, "fastp.{}".format(leftReadsBase))
                    newLeftReads.append(newLeft)
                newLeftReadsRepeats.append(",".join(newLeftReads))
            dataObj.ngsLeftReads = ";".join(newLeftReadsRepeats)
        elif dataObj.ngsRightReads and dataObj.ngsLeftReads == None:
            rightReadsRepeats = [i.strip() for i in dataObj.ngsRightReads.split(";")]
            newRightReadsRepeats = []
            for i in range(len(rightReadsRepeats)):
                rightReads = rightReadsRepeats[i].split(",")
                newRightReads = []
                for j in range(len(rightReads)):
                    rightReadsDir = os.path.dirname(rightReads[j])
                    rightReadsBase = os.path.basename(rightReads[j])
                    newRight = os.path.join(rightReadsDir, "fastp.{}".format(rightReadsBase))
                    newRightReads.append(newRight)
                newRightReadsRepeats.append(",".join(newRightReads))
            dataObj.ngsRightReads = ";".join(newRightReadsRepeats)
        else:
            raise Exception("The NGS data seem not to be single, please check it")

def validateFaAndFqFile(myFile):
    cmd = "seqkit head 100 {} | seqkit stat".format(myFile)
    validateOutput = os.popen(cmd)
    if "FASTQ" in validateOutput:
        return "fastq"
    elif "FASTA" in validateOutput:
        return "fasta"
    else:
        return False

def validateFile(myFile):
    if not os.path.exists(myFile):
        raise Exception("File '%s' not found! Please input again!" % myFile)

    if not os.path.isfile(myFile):
        raise Exception("File '%s' is not a file! Please input again!" % myFile)

    return True

def validateDir(myDir):
    if not os.path.exists(myDir):
        raise Exception("Dir '%s' not found! Please input again!" % myDir)

    if not os.path.isdir(myDir):
        raise Exception("Dir '%s' is not a directory! Please input again!" % myDir)

    return True