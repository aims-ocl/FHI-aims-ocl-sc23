#!/usr/bin/env python3

import os
import sys
import numpy
import scipy.lib.lapack.flapack as lapack

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
        




chf = open('hamiltonian.out','r')
cof = open('overlap-matrix.out','r')
phf = open('rrs_hamiltonian_center.out','r')
pof = open('rrs_overlap_center.out','r')

cs  = 82
ce  = 119
cd  = 200 
#cs  = 44
#ce  = 81
#cd  = 124

ps  = 1
pe  = 38
pd  = 38

cho = read_matrix(cd,chf,iop=0)
coo = read_matrix(cd,cof,iop=0)
for i in range(cd):
    for j in range(i):
        cho[i,j]=cho[j,i]
        coo[i,j]=coo[j,i]
ch  = cho[cs-1:ce,cs-1:ce]
co  = coo[cs-1:ce,cs-1:ce]
ph  = read_matrix(pd,phf,iop=1)[ps-1:pe,ps-1:pe]
po  = read_matrix(pd,pof,iop=1)[ps-1:pe,ps-1:pe]

dhami  = ch-ph
print('dhamt:%16.8f' %dhami.sum())
dovlp  = co-po
print('dovlp:%16.8f' %dovlp.sum())

thresh = 1.0E-5
print('     %5s%5s:%16s%16s' %('index','index','cluster','center'))
dhami  = 0.0
dovlp  = 0.0
for i in range(pd):
    for j in range(i):
        if abs(ch[i,j]-ph[i,j])>thresh:
            print('hamt:%5i%5i:%16.8f%16.8f' %(i+1,j+1,ch[i,j],ph[i,j]))
            dhami += abs(ch[i,j]-ph[i,j])
        if abs(co[i,j]-po[i,j])>thresh:
            print('ovlp:%5i%5i:%16.8f%16.8f' %(i+1,j+1,co[i,j],po[i,j]))
            dovlp += abs(co[i,j]-po[i,j])
print('dhamt:%16.8f' %dhami)
print('dovlp:%16.8f' %dovlp)


w,v,info = lapack.zhegv(ph,po,itype=1,compute_v=1,lower=0,lwork=2*pd-1,overwrite_a=0,overwrite_b=0)
print(w,'\n Diag. Info. ',info)


phff= open('rrs_hamiltonian_final.out','r')
poff= open('rrs_overlap_final.out','r')
phf = read_matrix(pd,phff,iop=1)[ps-1:pe,ps-1:pe]
pof = read_matrix(pd,poff,iop=1)[ps-1:pe,ps-1:pe]

#print_matrix(pd,phf)

w,v,info = lapack.zhegv(phf,pof,itype=1,compute_v=1,lower=0,lwork=2*pd-1,overwrite_a=0,overwrite_b=0)
print(w,'\n Diag. Info. ',info)

Gauphff= open('Gau_PBC_hamiltonian_final.out','r')
Gaupoff= open('Gau_PBC_overlap_final.out','r')
Gauphf = read_matrix(pd,Gauphff,iop=1)[ps-1:pe,ps-1:pe]
Gaupof = read_matrix(pd,Gaupoff,iop=1)[ps-1:pe,ps-1:pe]

#print_matrix(pd,phf)

w,v,info = lapack.zhegv(Gauphf,Gaupof,itype=1,compute_v=1,lower=0,lwork=2*pd-1,overwrite_a=0,overwrite_b=0)
print(w,'\n Diag. Info. ',info)

print('     %5s%5s:%16s%16s' %('index','index','cluster','center'))
dhami  = 0.0
dovlp  = 0.0
for i in range(pd):
    for j in range(i):
        if abs(Gauphf[i,j]-phf[i,j])>thresh:
            print('hamt:%5i%5i:%16.8f%16.8f' %(i+1,j+1,Gauphf[i,j],phf[i,j]))
            dhami += abs(Gauphf[i,j]-phf[i,j])
        if abs(Gaupof[i,j]-pof[i,j])>thresh:
            print('ovlp:%5i%5i:%16.8f%16.8f' %(i+1,j+1,Gaupof[i,j],pof[i,j]))
            dovlp += abs(Gaupof[i,j]-pof[i,j])
print('dhamt:%16.8f' %dhami)

dhami  = Gauphf-phf
print('dhamt:%16.8f' %dhami.sum())
dovlp  = Gaupof-pof
print('dovlp:%16.8f' %dovlp.sum())
