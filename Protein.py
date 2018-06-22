import numpy
import os
directoryname = "proteins!"
files = os.listdir(directoryname)
allproteins = set()

for filename in files:
   #TODO: os.path
   filepath = os.path.join(directoryname,filename)
   filestream = open(filepath, "r")
   for line in filestream:
       line = line.strip()
       allproteins.add(line)

allproteins = list(allproteins)

matrix = []
realsequence = []

for filename in files:
    filepath = os.path.join(directoryname,filename)
    sequence = open(filepath, "r")
    for line in sequence:
       line = line.strip()
       realsequence.append(line)
       vector = []
       for x in allproteins:
           cnt = realsequence.count(x)
           vector.append(cnt)
    matrix.append(vector)
    vector = []
    realsequence = []

matrix = numpy.array(matrix)
final = open("matrix.bin", "w")
matrix.tofile(final, "bw")
final.close()
