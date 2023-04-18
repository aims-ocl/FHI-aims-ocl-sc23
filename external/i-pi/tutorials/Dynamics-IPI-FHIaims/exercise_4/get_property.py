#!/usr/bin/python
import sys
import random

datatype = sys.argv[1]
name = sys.argv[2]
file0=open(name,'r')

# Pass the first 13 lines, which only indicate which information are printed
for i in range(13):
  file0.readline()

# Go through the whole file, line by line, and print the temperature against time
x1=1 # Second column - corresponds to the simulation time
if (datatype=="conserved"):
  x2=2
  output=open('conserved_'+name+'.dat', 'w')
elif (datatype=="temperature"):
  x2=3
  output=open('mean_temperature_'+name+'.dat', 'w')
elif (datatype=="kinetic"):
  x2=4
  output=open('kinetic_'+name+'.dat', 'w')
elif (datatype=="potential"):
  x2=5
  output=open('potential_'+name+'.dat', 'w')

while True:
  line=file0.readline()
  if not line:
    break
  #print line.split()[x1], line.split()[x2]
  #print line.split()[x2]
  output.write(line.split()[x2]+'\n')

output.close()

sentences=["Oh, you executed this script beautifully...", "You look beautiful today.", "Aren't you tired to launch scripts all day ?","I am still hungry, give me more inputs !","I am a rebel script. I refuse to be executed like this.","lalala lala la lalala"]
alea = random.randint(0,len(sentences)-1)
print sentences[alea]
