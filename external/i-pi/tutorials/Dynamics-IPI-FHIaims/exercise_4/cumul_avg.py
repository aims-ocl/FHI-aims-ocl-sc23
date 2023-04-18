#!/usr/bin/python
import sys

name = sys.argv[1]
file0=open(name,'r')

output=open('cumul_avg_'+name, 'w')

tot=0.0
counter=0
while True:
  line=file0.readline()
  if not line:
    break
  #print line.split()[x1], line.split()[x2]
  #print line.split()[x2]
  counter=counter+1
  tot=tot+float(line.split()[0])
  avg=tot/counter
  output.write(str(avg)+'\n')

output.close()

