#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: identify_as.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-05-12 22:07:52
Last modified: 2021-05-12 22:07:52
'''

def identify_as(dataObj=None, refParams=None, dirSpec=None):
    from find_as import find_as
    find_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
    from find_pa import find_pa
    find_pa(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
    from charaterize_as import charaterize_as
    charaterize_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)