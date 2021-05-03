#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: relationshipBetweenAlleleSpecificAndAS.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-05-14 20:55:41
Last modified: 2020-05-14 20:55:41
'''
import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from commonFuncs import *

def relationshipBetweenAlleleSpecifiAndAS(partial=False, nopartial=False, isoPairs=None):
    resultDict = {"partial": {}, "nopartial": {}}
    for isoIndex in range(len(isoPairs)):
        # partial
        if partial:
            partialHaplotype = "phased.partial.cleaned.human_readable.txt"
            if os.path.exists(partialHaplotype):
                partialDF = pd.read_csv(partialHaplotype, skiprows=[0], index_col=0, sep="\t")
                indexPairs = list(itertools.combinations(partialDF.index, 2))
                if len(set(isoPairs[isoIndex]) & set(partialDF.columns)) == 2:
                    for pairIndex in range(len(indexPairs)):
                        subPartialDF = partialDF.reindex(index=indexPairs[pairIndex], columns=isoPairs[isoIndex])
                        if (pd.DataFrame.max(subPartialDF) >= 5).all() and (pd.DataFrame.min(subPartialDF) <= 3).all():
                            partialChi2Results = chi2_contingency(subPartialDF)
                            partialChi2Pvalue = partialChi2Results[1]
                            if partialChi2Pvalue <= 0.001:
                                combination = "_".join(map(str, ["partial", isoIndex, pairIndex]))
                                rows = subPartialDF.index
                                columns = subPartialDF.columns
                                resultDict["partial"].update({combination: dict(zip(rows, columns))})


        # no-partial
        if nopartial:
            nopartialHaplotype = "phased.nopartial.cleaned.human_readable.txt"
            if os.path.exists(nopartialHaplotype):
                nopartialDF = pd.read_csv(nopartialHaplotype, skiprows=[0], index_col=0, sep="\t")
                indexPairs = list(itertools.combinations(nopartialDF.index, 2))
                if len(set(isoPairs[isoIndex]) & set(nopartialDF.columns)) == 2:
                    for pairIndex in range(len(indexPairs)):
                        subNopartialDF = nopartialDF.reindex(index=indexPairs[pairIndex], columns=isoPairs[isoIndex])
                        nopartialChi2Results = chi2_contingency(subNopartialDF)
                        nopartialChi2Pvalue = nopartialChi2Results[1]
                        if nopartialChi2Pvalue <= 0.05:
                            combination = "_".join(map(str, ["nopartial", isoIndex, pairIndex]))
                            rows = subNopartialDF.index
                            columns = subNopartialDF.columns
                            resultDict["nopartial"].update({combination: dict(zip(rows, columns))})

        if not partial and not nopartial:
            raise Exception("You must specify either partial or nopartial file")
    return resultDict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-i", type=str,
                        help="All isoforms from iFLAS pipeline")
    parser.add_argument("-f", type=str,
                        help="The flnc reads in bed format from iFLAS pipeline")
    parser.add_argument("-a", type=str,
                        help="AS file from iFLAS pipeline")
    parser.add_argument("-c", type=str,
                        help="The collapsed isoform file from cDNA_Cupcake")
    parser.add_argument("-o", type=str,
                        help="Output file")
    parser.add_argument("-asType", type=str,
                        help="The alternative splicing type")
    parser.add_argument("-filterByCount", type=int, default=0,
                        help="Filter the collapsed isoforms by the count of supporting flnc reads")
    args = parser.parse_args()
    relationshipBetweenAlleleSpecifiAndAS()
