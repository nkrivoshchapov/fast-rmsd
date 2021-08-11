import sys,os
from copy import deepcopy

refdata = []
reflines = open("results_ref.dat", "r").readlines()
for line in reflines:
    curdata = []
    parts = line.replace("\n", "").split()
    for part in parts:
        curdata.append(float(part))
    refdata.append(deepcopy(curdata))

cppdata = []
cpplines = open("results_cpp.dat", "r").readlines()
for line in cpplines:
    curdata = []
    parts = line.replace("\n", "").split()
    for part in parts:
        curdata.append(float(part))
    cppdata.append(deepcopy(curdata))

for i in range(1000):
    for j in range(1000):
        if abs(cppdata[i][j] - refdata[i][j]) > 0.0001:
            print("%d & %d" % (i, j))
        # else:
        #     print("OK : %d & %d" % (i, j))
