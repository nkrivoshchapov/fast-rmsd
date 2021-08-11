from __future__ import print_function
import numpy as np
import rmsd, ntpath, glob
from copy import deepcopy
import progressbar

logger = rmsd.createLogger("Main")
def getxyz(filename):
    xyz = []
    lines = open(filename,"r").readlines()
    for i in range(2,len(lines)):
        vector = np.array([float(lines[i].split()[1]),
                           float(lines[i].split()[2]),
                           float(lines[i].split()[3])])
        xyz.append(vector)
    return np.array(xyz)

mols = []
filelist = glob.glob("../testset/*xyz")
filelist.sort()
for file in filelist:
    print("Reading " + ntpath.basename(file))
    mols.append(getxyz(file))
    mols[len(mols)-1] -= rmsd.centroid(mols[len(mols)-1])

rmsds = []
bar = progressbar.ProgressBar(max_value=len(filelist))
for i in range(len(filelist)):
    curline = []
    for j in range(len(filelist)):
        U_full = rmsd.kabsch(mols[i], mols[j])
        Temp_full = np.dot(mols[i], U_full)
        curline.append(str(rmsd.rmsd(Temp_full, mols[j])))
    rmsds.append(" ".join(curline))
    bar.update(i)

wfile = open("../results_ref.dat", "w")
wfile.write("\n".join(rmsds))
wfile.close()

# A = getxyz("molC1.xyz")
# logger.info("A = \n%s" % repr(A))
# logger.info("centroid of A = \n%s" % repr(rmsd.centroid(A)))
# A -= rmsd.centroid(A)
# logger.info("A - centroid(A) = \n%s" % repr(A))

# B = getxyz("molD1.xyz")
# logger.info("B = \n%s" % repr(B))
# logger.info("centroid of B = \n%s" % repr(rmsd.centroid(B)))
# B -= rmsd.centroid(B)
# logger.info("B - centroid(B) = \n%s" % repr(B))

# U_full = rmsd.kabsch(A, B)
# logger.info("kabsch(A,B) = \n%s" % repr(U_full))
# Temp_full = np.dot(A, U_full)
# logger.info("np.dot(A, kabsch(A,B)) = \n%s" % repr(Temp_full))
# logger.info("RMSD = " + repr(rmsd.rmsd(Temp_full, B)))
