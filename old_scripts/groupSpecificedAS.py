#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: groupSpecificedAS.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-06-25 23:34:30
Last modified: 2020-06-25 23:34:30
'''

asTypes = ["IR", "SE", "A3SS", "A5SS"]

def processAsFile(asFile):
    myDict = {}
    with open(asFile) as f:
        for line in f.readlines():
            myDict[line.strip().split("\t")[3]] = line.strip()
    return myDict

def main():
    group1dirs = ["", ""]
    group2dirs = ["", ""]
    for i in asTypes:
        out1 = open("group1specific.{}.lst".format(i), "w")
        out2 = open("group2specific.{}.lst".format(i), "w")
        if i == "SE":
            extension = "bed12+"
        else:
            extension = "bed6+"
        group1File1 = "{}/{}.{}".format(group1dirs[0], i.lower(), extension)
        group1File2 = "{}/{}.{}".format(group1dirs[1], i.lower(), extension)
        group2File1 = "{}/{}.{}".format(group2dirs[0], i.lower(), extension)
        group2File2 = "{}/{}.{}".format(group2dirs[1], i.lower(), extension)

        group1File1_as = processAsFile(group1File1)
        group1File2_as = processAsFile(group1File2)
        group2File1_as = processAsFile(group2File1)
        group2File2_as = processAsFile(group2File2)

        group1intersect = set(group1File1_as.keys()) & set(group1File2_as.keys())
        group2union = set(group2File1_as.keys()) | set(group2File2_as.keys())

        group1union = set(group1File1_as.keys()) | set(group1File2_as.keys())
        group2intersect = set(group2File1_as.keys()) & set(group2File2_as.keys())

        for j in group1intersect:
            if j not in group2union:
                print >> out1, j

        for j in group2intersect:
            if j not in group1union:
                print >> out2, j
        out1.close()
        out2.close()

if __name__ == '__main__':
    main()