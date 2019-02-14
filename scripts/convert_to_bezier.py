import numpy as np
import sys
import csv
import math
from numpy.linalg import inv
import argparse

parser = argparse.ArgumentParser(description="convert traj to bezier")
parser.add_argument("inputfile", help="input file destination")
parser.add_argument("outputfile", help="output file destination")
args = parser.parse_args()

trajfile = args.inputfile
outputfile = args.outputfile


def comb(a, b):
    f = math.factorial
    return f(a) / f(b) / f(a-b)


trajstr = []

with open(trajfile, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        trajstr.append(row)

trajstr = trajstr[1:]
traj = []
for l in trajstr:
    tr = []
    for p in l[:25]:
        tr.append(float(p))
    traj.append(tr)


matrix = []
for k in range(0, 8):
    matrixrow = []
    for i in range(0, 8):
        if i > k:
            matrixrow.append(0)
        else:
            eff = comb(7-i, k-i) * pow(-1, k-i) * comb(7,i)
            matrixrow.append(eff)
    matrix.append(matrixrow)



matrixnp = np.array(matrix)
invmtr = inv(matrixnp)

out = open(outputfile, "w")

out.write("#duration,p0x,p0y,p0z,...\n")

for trj in traj:
    duration = trj[0]
    multip = 1
    for i in range(8):
        trj[1+i] *= multip
        trj[9+i] *= multip
        trj[17+i] *= multip
        multip *= duration

    xnp = np.transpose(np.array(trj[1:9]))
    ynp = np.transpose(np.array(trj[9:17]))
    znp = np.transpose(np.array(trj[17:]))
    xctrl = np.matmul(invmtr, xnp).tolist()
    yctrl = np.matmul(invmtr, ynp).tolist()
    zctrl = np.matmul(invmtr, znp).tolist()
    pstr = ""
    pstr += str(duration) + ","
    for i in range(7):
        pstr += str(xctrl[i]) + ","
        pstr += str(yctrl[i]) + ","
        pstr += str(zctrl[i]) + ","
    pstr += str(xctrl[-1]) + ","
    pstr += str(yctrl[-1]) + ","
    pstr += str(zctrl[-1])
    out.write(pstr + "\n")
out.close()
