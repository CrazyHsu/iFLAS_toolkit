#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: gtfObjs.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-07-27 14:53:18
Last modified: 2019-07-27 14:53:19
'''
from GeneModel import *

PARENT, ID, NAME = "parent", "id", "name"
TRANS_TYPE, EXON_TYPE, CDS_TYPE, UTR_TYPE, SELENO_TYPE = "transcript", "exon", "cds", "utr", "selenocysteine"
START_CODON, STOP_CODON = "start_codon", "stop_codon"

F_UTR, T_UTR = "five_prime_utr", "three_prime_utr"

GENE_ID_ATTR, TRANS_ID_ATTR, EXON_ID_ATTR, PROTEIN_ID_ATTR = \
    "gene_id", "transcript_id", "exon_id", "protein_id"
EXON_NO_ATTR, GENE_NAME_ATTR, TRANS_NAME_ATTR = \
    "exon_number", "gene_name", "transcript_name"
GENE_BIOTYPE, TRANS_BIOTYPE = "gene_biotype", "transcript_biotype"
GENE_SOURCE, TRANS_SOURCE = "gene_source", "transcript_source"

ID_ATTR = [GENE_ID_ATTR, TRANS_ID_ATTR]
NAME_ATTR = [GENE_NAME_ATTR, TRANS_NAME_ATTR]
GENE_OPT_ATTR = [GENE_NAME_ATTR, GENE_BIOTYPE, GENE_SOURCE]
TRANS_OPT_ATTR = [TRANS_NAME_ATTR, TRANS_BIOTYPE, TRANS_SOURCE]
NAME2ID = {GENE_NAME_ATTR: GENE_ID_ATTR, TRANS_NAME_ATTR: TRANS_ID_ATTR}

class GTF_Line(object):
    def __init__(self, line):
        if line.startswith("#"):
            raise ValueError()
        self.fields = line.strip().split("\t")
        if len(self.fields) != 9:
            raise ValueError("This is incorrect GTF format where the correct GTF format should contain 9 column")
        self.setFieldInfo()
        self.attrs = self.setAttrs()
        self.setAttrInfo()

    def setAttrs(self):
        attrList = self.fields[8].split(";")
        attrDict = {}
        for p in attrList:
            if not p: continue
            pair = p.split()
            key = pair[0].strip()
            value = pair[1].strip().replace('"', '')
            attrDict[key] = value.upper()

        absentList = []
        for attr in ID_ATTR:
            if attr not in attrDict:
                if self.feature.upper() == "GENE":
                    attrDict["transcript_id"] = None
                else:
                    absentList.append(attr)
        if len(absentList) > 0:
            raise ValueError("Please check if you have input the right GTF file with {} attribute".format(", ".join(absentList)))

        for attr in NAME_ATTR:
            if attr not in attrDict:
                attrDict[attr] = attrDict[NAME2ID[attr]]
        return attrDict

    def setFieldInfo(self):
        self.chrId = self.fields[0]
        self.source = self.fields[1]
        self.feature = self.fields[2]
        self.start = int(self.fields[3])
        self.end = int(self.fields[4])
        self.score = self.fields[5]
        self.strand = self.fields[6]
        self.phase = self.fields[7]

    def setAttrInfo(self):
        self.geneId = self.attrs[GENE_ID_ATTR]
        self.transId = self.attrs[TRANS_ID_ATTR]

    def __str__(self):
        return "%s:%d-%d,%s" % (self.chrId, self.start, self.end, self.feature)

class GTFLoader(object):
    def __init__(self, gtfTabix, **args):
        # self.model = self.loadGeneModels(path, **args)
        model = self.loadGeneModels(gtfTabix, **args)
        self.allGenes = model.allGenes
        self.allChrom = model.allChr
        self.geneModel = model.geneModel

    def loadGeneModels(self, gtfTabix, **args):
        geneModel = GeneModel()
        for line in gtfTabix.fetch():
        # with open(path) as f:
        #     for line in f.readlines():
            try:
                record = GTF_Line(line.strip("\n"))
            except:
                continue
            record.chrId = record.chrId.lower()
            chrId = record.chrId
            feature = record.feature
            strand = record.strand
            start = record.start
            end = record.end
            geneModel.updateChrom(1, end, chrId)
            gene = self.updateGene2Model(geneModel, record, start, end)

            transId = record.attrs[TRANS_ID_ATTR]
            feat = feature.lower()
            if feat == TRANS_TYPE:
                self.updateTrans2Gene(gene, record)
            elif feat == EXON_TYPE:
                self.updateExon2Gene(gene, record)
            elif feat == CDS_TYPE:
                self.updateCDS2Gene(gene, record)
            elif UTR_TYPE in feat:
                self.updateUTR2Gene(gene, record)
            elif feat == START_CODON:
                gene.startCodons[transId] = (start, end)
                if transId in gene.trans:
                    gene.trans[transId].startCodon = (start, end)
            elif feat == STOP_CODON:
                gene.stopCodons[transId] = (start, end)
                if transId in gene.trans:
                    gene.trans[transId].stopCodon = (start, end)
        self.fillGeneModel(geneModel)
        return geneModel

    def fillGeneModel(self, geneModel):
        for geneId, geneObj in geneModel.allGenes.items():
            if len(geneObj.utrs) == 0 and len(geneObj.cds) == 0:
                continue
            elif len(geneObj.utrs) == 0 and len(geneObj.cds) != 0:
                for transId, transObj in geneObj.trans.items():
                    # futrStart = transObj.minpos if transObj.strand == "+" else transObj.startCodon[1] + 1
                    # tutrEnd = transObj.maxpos if transObj.strand == "+" else transObj.stopCodon[1]
                    # futrEnd = transObj.startCodon[0] - 1 if transObj.strand == "+" else transObj.maxpos
                    # tutrStart = transObj.stopCodon[0] if transObj.strand == "+" else transObj.minpos

                    if transObj.startCodon:
                        futrStart = transObj.minpos if transObj.strand == "+" else transObj.startCodon[1] + 1
                        futrEnd = transObj.startCodon[0] - 1 if transObj.strand == "+" else transObj.maxpos
                    else:
                        futrStart = transObj.minpos if transObj.strand == "+" else max([i.end for i in transObj.cds]) + 1
                        futrEnd = min([i.start for i in transObj.cds]) - 1 if transObj.strand == "+" else transObj.maxpos
                    if transObj.stopCodon:
                        tutrStart = transObj.stopCodon[0] if transObj.strand == "+" else transObj.minpos
                        tutrEnd = transObj.maxpos if transObj.strand == "+" else transObj.stopCodon[1]
                    else:
                        tutrStart = max([i.end for i in transObj.cds]) + 1 if transObj.strand == "+" else transObj.minpos
                        tutrEnd = transObj.maxpos if transObj.strand == "+" else min([i.start for i in transObj.cds]) - 1
                    # fcds = transObj.cds[0] if transObj.strand == "+" else transObj.cds[-1]
                    # tcds = transObj.cds[-1] if transObj.strand == "+" else transObj.cds[0]

                    futrs, tutrs = [], []
                    if transObj.strand == "+":
                        for exon in transObj.exons:
                            futr, tutr = None, None
                            if exon.start >= futrStart and exon.end <= futrEnd:
                                futr = UTR(exon.start, exon.end, transObj.chrId, transObj.strand, F_UTR, attr={})
                            elif exon.start <= futrEnd and exon.end > futrEnd:
                                futr = UTR(exon.start, futrEnd, transObj.chrId, transObj.strand, F_UTR, attr={"fiveJunc": 1})
                            if exon.start >= tutrStart and exon.end <= tutrEnd:
                                tutr = UTR(exon.start, exon.end, transObj.chrId, transObj.strand, T_UTR, attr={})
                            elif exon.start < tutrStart and exon.end > tutrStart:
                                tutr = UTR(tutrStart, exon.end, transObj.chrId, transObj.strand, T_UTR, attr={"threeJunc": 1})
                            if futr != None: futrs.append(futr)
                            if tutr != None: tutrs.append(tutr)
                    else:
                        for exon in transObj.exons:
                            futr, tutr = None, None
                            if exon.start >= futrStart and exon.end <= futrEnd:
                                futr = UTR(exon.start, exon.end, transObj.chrId, transObj.strand, F_UTR, attr={})
                            elif exon.start <= futrStart and exon.end > futrStart:
                                futr = UTR(futrStart, exon.end, transObj.chrId, transObj.strand, F_UTR, attr={"fiveJunc": 1})
                            elif exon.start >= tutrStart and exon.end <= tutrEnd:
                                tutr = UTR(exon.start, exon.end, transObj.chrId, transObj.strand, T_UTR, attr={})
                            elif exon.start < tutrEnd and exon.end >=tutrEnd:
                                tutr = UTR(exon.start, tutrEnd, transObj.chrId, transObj.strand, T_UTR, attr={"threeJunc": 1})
                            if futr != None: futrs.append(futr)
                            if tutr != None: tutrs.append(tutr)
                    for utr in futrs + tutrs:
                        if (utr.start, utr.end) not in geneObj.utrDict:
                        #     u = geneObj.utrDict[(utr.start, utr.end)]
                        #     u.feature = utr.feature
                        # else:
                            geneObj.utrs.append(utr)
                            geneObj.utrDict[(utr.start, utr.end)] = utr
                        geneObj.trans[transId].updateUTR(utr)

                    # futr = UTR(futrStart, futrEnd, transObj.chrId, transObj.strand, "five_prime_utr", attr={})
                    # tutr = UTR(tutrStart, tutrEnd, transObj.chrId, transObj.strand, "three_prime_utr", attr={})
                    # for utr in [futr, tutr]:
                    #     if (utr.start, utr.end) in geneObj.utrDict:
                    #         u = geneObj.utrDict[(utr.start, utr.end)]
                    #         u.feature = utr.feature
                    #     else:
                    #         geneObj.utrs.append(utr)
                    #         geneObj.utrDict[(utr.start, utr.end)] = utr
                    #         geneObj.trans[transObj.transId].updateUTR(utr)

    def updateCDS2Gene(self, gene, record):
        start = record.start
        end = record.end
        chrId = record.chrId
        strand = record.strand
        phase = record.phase
        cdsObj = CDS(start, end, chrId, strand, phase, record.attrs)
        cdsBound = (record.start, record.end)
        if cdsBound in gene.cdsDict:
            cds = gene.cdsDict[cdsBound]
        else:
            cds = cdsObj
            gene.cds.append(cds)
            gene.cdsDict[cdsBound] = cds
        gene.cds = sorted(gene.cds, key=lambda n: n.start)
        trans = self.updateTrans2Gene(gene, record)
        trans.updateCDS(cdsObj)

    def updateGene2Model(self, geneModel, record, start, end):
        minpos = min(start, end)
        maxpos = max(start, end)
        chrId = record.chrId
        geneId = record.attrs[GENE_ID_ATTR]
        strand = record.strand
        try:
            gene = geneModel.geneModel[chrId][geneId]
            gene.minpos = min(minpos, gene.minpos)
            gene.maxpos = max(maxpos, gene.maxpos)
        except:
            geneOptAttrs = dict([(n, record.attrs[n]) for n in record.attrs if n in GENE_OPT_ATTR])
            gene = Gene(geneId, start, end, chrId, strand, record.attrs[GENE_NAME_ATTR], geneOptAttrs)
            geneModel.geneModel[chrId][geneId] = gene
            geneModel.allGenes[geneId] = gene
        return gene

    def updateExon2Gene(self, gene, record):
        exonObj = Exon()

        chrId = record.chrId
        end = record.end
        strand = record.strand
        start = record.start
        transId = record.attrs[TRANS_ID_ATTR]
        transName = record.attrs[TRANS_NAME_ATTR]
        exonAttr = record.attrs

        exonObj.start = start
        exonObj.end = end
        exonObj.strand = strand
        exonObj.exonAttr = exonAttr
        exonObj.chrId = chrId

        if transId in gene.trans:
            newTrans = gene.trans[transId]
        else:
            transAttr = {PARENT: gene.geneId, ID: transId, NAME: transName}
            newTrans = Trans(transId, start, end, chrId, strand, transAttr)

        exonBound = (start, end)
        if exonBound in gene.exonDict:
            exon = gene.exonDict[exonBound]
        else:
            exon = exonObj
            gene.exons.append(exon)
            gene.exonDict[exonBound] = exon

        gene.exons = sorted(gene.exons, key=lambda n: n.start)
        trans = gene.updateTrans(newTrans)
        trans.updateExons(exon, transId=transId)
        # trans.updateExons(exon)
        # return exon

    def updateSeleno2Gene(self, gene, record):
        chrId = record.chrId
        strand = record.strand
        start = record.start
        end = record.end
        transId = record.attrs[TRANS_ID_ATTR]

    def updateTrans2Gene(self, gene, record):
        chrId = record.chrId
        strand = record.strand
        start = record.start
        end = record.end
        transId = record.attrs[TRANS_ID_ATTR]
        transName = record.attrs[TRANS_NAME_ATTR]
        try:
            return gene.trans[transId]
        except:
            transAttr = {PARENT: gene.geneId, ID: transId, NAME: transName}
            gene.updateTrans(Trans(transId, start, end, chrId, strand, transAttr))
            return gene.trans[transName]

    def updateUTR2Gene(self, gene, record):
        chrId = record.chrId
        start = record.start
        end = record.end
        strand = record.strand
        transId = record.attrs[TRANS_ID_ATTR]
        transName = record.attrs[TRANS_NAME_ATTR]
        utrType = record.feature
        utrAttr = record.attrs
        utrObj = UTR(start, end, chrId, strand, utrType, utrAttr)
        utrBound = (start, end)
        tmp = self.updateTrans2Gene(gene, record)
        if utrBound in gene.utrDict:
            utr = gene.utrDict[utrBound]
            utr.feature = utrType
        else:
            gene.utrs.append(utrObj)
            gene.utrDict[utrBound] = utrObj
            gene.trans[transId].updateUTR(utrObj)
        gene.utrs = sorted(gene.utrs, key=lambda n: n.start)

        # gene.updateUTR(UTR(start, end, chrId, strand, utrType, utrAttr))
        # startCodon = gene.startCodons[transId] if transId in gene.startCodons else None
        # endCodon = gene.endCodons[transId] if transId in gene.endCodons else None
