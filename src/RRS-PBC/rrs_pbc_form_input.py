#!/usr/bin/env python3

import os
import sys
import numpy
import scipy

def read_matrix(dimension, flow, iop=0) :
    mat = numpy.matrix(numpy.arange(dimension**2).reshape((dimension,dimension)),dtype=numpy.float128)
    m = numpy.mod(dimension,5)
    n = int((dimension-m)/5)
    for i in range(n):
        tmpline = flow.readline()
        if i!=0 and iop==0:
            tmpline = flow.readline()
        for j in range(dimension):
            if iop==0:
                tmpline = flow.readline().strip().split()[1:]
            else:
                tmpline = flow.readline().strip().split()[2:]
            ms          = i*5
            me          = (i+1)*5
            mat[j,ms:me]=[float(tmpline[k]) for k in range(5)]
    if m!=0:
        tmpline = flow.readline()
        if iop==0:
            tmpline = flow.readline()
        for j in range(dimension):
            if iop==0:
                tmpline = flow.readline().strip().split()[1:]
            else:
                tmpline = flow.readline().strip().split()[2:]
            ms          = dimension-m
            me          = dimension
            mat[j,ms:me]=[float(tmpline[k]) for k in range(m)]
    return mat

def print_matrix(dimension, mat) :
    m = numpy.mod(dimension,5)
    n = int((dimension-m)/5)
    for i in range(n):
        print('     '+'%16i'*5 %(tuple(range(i*5+1,(i+1)*5+1,1))))
        for j in range(dimension):
            ms          = i*5
            me          = (i+1)*5
            print('%5i' %(j+1)+'%16.8f'*5 %(tuple(mat.getA()[j,ms:me])))
    if m!=0:
        print('     '+'%16i'*m %(tuple(range(n*5+1,dimension+1,1))))
        for j in range(dimension):
            ms          = dimension-m
            me          = dimension
            print('%5i' %(j+1)+'%16.8f'*m %(tuple(mat.getA()[j,ms:me])))
    return
        
cf = open('control.in','r')
cs = cf.readlines()
cf.close()
gf = open('geometry.in','r')
gs = gf.readlines()
gf.close()

# loading geometry information
geometry = []
for x in gs:
    tmplist = x.strip().


