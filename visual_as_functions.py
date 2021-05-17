from commonFuncs import *
# from commonObjs import *
import cPickle as pickle
import numpy as np

class MainSection(object):
    '''
    Set values of variables in main section.
    '''
    def __init__(self, **args):
        self.sectionName = getAttribute("sectionName", "[SinglePlotConfig]", **args)
        self.legend      = getAttribute("legend", True, **args)
        self.width       = getAttribute("width", 9.0, **args)
        self.height      = getAttribute("height", 11.0, **args)
        self.fontsize    = getAttribute("fontsize", 12, **args)
        self.output_file = getAttribute("fout", None, **args)
        self.shrink_introns = getAttribute("shrink_introns", False, **args)

    def printStr(self):
        return "%s\nlegend = %s\nwidth = %.1f\nheight = %.1f\noutput_file = %s\nshrink_introns = %s\n" % \
               (self.sectionName, self.legend, self.width, self.height, self.output_file, self.shrink_introns)

class PlotSection(object):
    '''
    Set values of variables in plot section.
    '''
    def __init__(self, **args):
        self.section_name = getAttribute("section_name", "[testGeneModel]", **args)
        self.plot_type = getAttribute("plot_type", "gene", **args)
        self.relative_size = getAttribute("relative_size", 10.0, **args)
        self.source_file = getAttribute("source_file", None, **args)
        self.title_string = getAttribute("title_string", "testTitle", **args)
        if self.plot_type == "gene":
            self.gene_name = getAttribute("gene_name", "testGeneName", **args)
            self.file_format = getAttribute("file_format", "gene_model", **args)
            self.hide = getAttribute("hide", False, **args)

    def printStr(self):
        returnStr = "%s\nplot_type = %s\nrelative_size = %.1f\nsource_file = %s\ntitle_string = %s\n" % \
                    (self.section_name, self.plot_type, self.relative_size, self.source_file, self.title_string)
        if self.plot_type == "gene":
            returnStr += "gene_name = %s\nfile_format = %s\nhide = %s\n" % (self.gene_name, self.file_format, self.hide)
        return returnStr
    
#####################################
def getAttribute(key, default, **args):
    return default if key not in args else args[key]

def getFileRowCounts(inFile):
    return len(open(inFile).readlines())

def resizeTrackRatio(itemCounts):
    originRelativeSize = 10.0
    if itemCounts <= 10:
        return originRelativeSize
    elif 10 <= itemCounts and itemCounts <= 100:
        return 2 * originRelativeSize * np.log10(itemCounts) - originRelativeSize
    else:
        return originRelativeSize * 3

def processTargetGenes(targetGenes):
    if os.path.isfile(targetGenes):
        return getDictFromFile(targetGenes)
    else:
        return targetGenes.split(",")
    # elif "," in targetGenes:
    #     return targetGenes.split(",")
    # else:
    #     raise Exception("Please input target genes in corrected format: "
    #                     "1. Input the genes in string split by comma; "
    #                     "2. Input the genes in a file one gene per line")

def parallelPlotterAnno(gene, gpeTargetGenePickle, sampleTargetGenePickle, dataObjs, dirSpec):
    projectName, sampleName = dataObjs.project_name, dataObjs.sample_name
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    isoViewDir = os.path.join(baseDir, "isoViewer")
    resolveDir(os.path.join(isoViewDir, gene))

    sampleTargetGeneObj = pickle.loads(sampleTargetGenePickle)
    targetGeneChrom = sampleTargetGeneObj.chrom
    targetGeneMinpos = sampleTargetGeneObj.minpos
    targetGeneMaxpos = sampleTargetGeneObj.maxpos
    targetGeneStrand = sampleTargetGeneObj.strand

    # Gene model section
    geneModelGPE = "{}.gpe".format(gene)
    geneModelGTF = "{}.gtf".format(gene)
    geneModelGFF = "{}.gff".format(gene)

    gpeTargetGene = pickle.loads(gpeTargetGenePickle)
    geneModelGPEOut = open(geneModelGPE, "w")
    gpeMinposList, gpeMaxposList = [targetGeneMinpos], [targetGeneMaxpos]
    for i in gpeTargetGene:
        gpeMinposList.append(i.txStart)
        gpeMaxposList.append(i.txEnd)
        print >> geneModelGPEOut, i
    geneModelGPEOut.close()
    plotMinpos = min(gpeMinposList)
    plotMaxpos = max(gpeMaxposList)
    targetGeneRegion = "{}:{}-{}".format(targetGeneChrom, plotMinpos, plotMaxpos)
    os.system("genePredToGtf file {} {} -source=iFLAS".format(geneModelGPE, geneModelGTF))
    os.system("gene_model_to_splicegraph.py -m {} -o {} 2>/dev/null".format(geneModelGTF, geneModelGFF))
    geneModelHidePlot = PlotSection(section_name="[GeneModelGraph]", source_file=geneModelGTF,
                                    gene_name=gene, relative_size=5.0,
                                    title_string="Gene Model for %gene", hide=True)
    geneModelItemCount = getFileRowCounts(geneModelGPE)
    geneModelRelativeSize = resizeTrackRatio(geneModelItemCount)
    geneModelVisiblePlot = PlotSection(section_name="[GeneModelIsoformsGraph]", plot_type="isoforms",
                                       source_file=geneModelGFF, relative_size=geneModelRelativeSize,
                                       title_string="Gene Model for %gene [{}({})]".format(
                                           targetGeneRegion, targetGeneStrand))

    # Corrected tgs reads from the pipeline
    postCorrIsoforms = sampleTargetGeneObj.reads
    postCorrIsoGPE = "{}.iFLAS.gpe".format(gene)
    postCorrIsoGTF = "{}.iFLAS.gtf".format(gene)
    postCorrIsoGFF = "{}.iFLAS.gff".format(gene)
    postCorrIsoGPEOut = open(postCorrIsoGPE, "w")
    for read in postCorrIsoforms:
        print >> postCorrIsoGPEOut, postCorrIsoforms[read].go_to_gpe()
    postCorrIsoGPEOut.close()
    os.system("genePredToGtf file {} {} -source=iFLAS".format(postCorrIsoGPE, postCorrIsoGTF))
    postCorrIsoItemCounts = len(postCorrIsoforms)
    postCorrIsoRelativeSize = resizeTrackRatio(postCorrIsoItemCounts)
    os.system("gene_model_to_splicegraph.py -m {} -o {} -a 2>/dev/null".format(postCorrIsoGTF, postCorrIsoGFF))
    postCorrIsoPlotType = "isoforms"
    postCorrIsoPlot = PlotSection(section_name="[AllReadsCollapse]", plot_type=postCorrIsoPlotType,
                                  source_file=postCorrIsoGFF, relative_size=postCorrIsoRelativeSize,
                                  title_string="Corrected isoforms and AS events in %s from TGS data" % gene)


    figOut = gene + ".pdf"
    cfgOut = open(gene + ".cfg", "w")
    majorItemCount = geneModelItemCount + postCorrIsoItemCounts
    figHeight = 10 if majorItemCount <= 50 else 15 if majorItemCount <= 150 else 20
    mainSec = MainSection(fout=figOut, height=figHeight)
    print >> cfgOut, mainSec.printStr()
    print >> cfgOut, geneModelHidePlot.printStr()
    print >> cfgOut, geneModelVisiblePlot.printStr()
    print >> cfgOut, postCorrIsoPlot.printStr()

    # Abundance from sam file in NGS pipeline
    if dataObjs.ngs_left_reads != None or dataObjs.ngs_right_reads != None:
        ngsSams = []
        for n in range(len(dataObjs.ngs_left_reads.split(";"))):
            repeatName = "repeat" + str(n)
            bamFile = os.path.join(baseDir, "mapping", "rna-seq", "alignment", "{}/{}.sorted.bam".format(repeatName, repeatName))
            targetSam = "{}.{}.sam".format(repeatName, gene)
            cmd = "samtools view -h {} {} > {}".format(bamFile, targetGeneRegion, targetSam)
            subprocess.call(cmd, shell=True)
            ngsSams.append(targetSam)
            readsPlot = PlotSection(section_name="[Reads_%s]" % (repeatName), plot_type="read_depth",
                                    source_file=targetSam, relative_size=5.0,
                                    title_string="%s Read Coverage in sample %s" % (gene, repeatName))
            print >> cfgOut, readsPlot.printStr()

    cfgOut.close()
    os.system("plotter.py %s.cfg 2>/dev/null" % gene)
    # removeFiles(ngsSams)
