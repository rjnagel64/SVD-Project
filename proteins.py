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
           if x in realsequence:
               vector.append("1")
           else:
               vector.append("0")
    matrix.append(vector)
    vector = []
    realsequence = []


print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in matrix]))