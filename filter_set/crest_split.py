import sys,os,glob
from copy import copy
for file in glob.glob("./7UPJ*crestdone.xyz"):
    lines = open(file,"r").readlines()
    start = lines[0]
    myname = file.split("/")[1].replace("_crestdone", "").replace(".xyz","")
    count = 0
    for i in range(len(lines)):
        if lines[i] == lines[0]:
            if i != 0:
                count += 1
                wfile = open("%s_%d.xyz" % (myname, count), "w")
                wfile.write("".join(sublines))
                wfile.close()
            sublines = []
        sublines.append(copy(lines[i]))
