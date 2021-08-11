import glob
from shutil import move

count = 0
for file in glob.glob("./*xyz"):
    move(file, ("%3d.xyz" % count).replace(" ","0"))
    count += 1
