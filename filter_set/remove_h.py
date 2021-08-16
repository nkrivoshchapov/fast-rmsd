import sys,os,glob

for file in glob.glob("./*xyz"):
    lines = open(file,"r").readlines()
    todel = []
    for i in range(2,len(lines)):
        if "H" in lines[i]:
            todel.append(i)
    count = int(lines[0])-len(todel)
    lines[0] = "%d\n" % count
    for i in reversed(todel):
        del lines[i]
    
    wfile = open(file.replace(".xyz","_noh.xyz"),"w")
    wfile.write("".join(lines))
    wfile.close()
