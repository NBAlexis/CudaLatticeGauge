import numpy as np
from matplotlib import pyplot as plt

from Distance.UsefulFunctions import findDistances
from Visualization import errorbar, MarkerStyle

allmarkers = [MarkerStyle.circle, MarkerStyle.triangle_up, MarkerStyle.plus, MarkerStyle.x]
folderstart = "G:\\"
for k in range(4):
    k1resn = np.loadtxt(f'{folderstart}\\WilsonLoopRes\\BS01DistanceNK{k+1}.csv', delimiter=',')
    k1resr = np.loadtxt(f'{folderstart}\\WilsonLoopRes\\BS01DistanceRK{k+1}.csv', delimiter=',')

    d = k1resn[0,:]
    resv1 = k1resn[1,:]
    ress1 = k1resn[2,:]

    resv2 = k1resn[41,:]
    ress2 = k1resn[42,:]

    resv3 = k1resr[1,:]
    ress3 = k1resr[2,:]

    resv4 = k1resr[41,:]
    ress4 = k1resr[42,:]

    errorbar(d, [resv1, resv2, resv3, resv4], [ress1, ress2, ress3, ress4],
             xlabel="$r/a$",
             ylabel="$\\langle W\\rangle$",
             markers=allmarkers,
             lineWidth=0.5,
             legends=["$\\langle W_{s.o.}\\rangle$ at $\\beta=5.3$","$\\langle W_{s.o.}\\rangle$ at $\\beta=5.8$","$\\langle W_{o.o.}\\rangle$ at $\\beta=5.3$","$\\langle W_{o.o.}\\rangle$ at $\\beta=5.3$"],
             savefile=f"../data/distance/BS01K{k+1}.pdf")

for k in range(4):
    k1resn = np.loadtxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceNK{k+1}.csv', delimiter=',')
    k1resr = np.loadtxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceRK{k+1}.csv', delimiter=',')

    d = k1resn[0,:]
    resv1 = k1resn[1,:]
    ress1 = k1resn[2,:]

    resv2 = k1resn[41,:]
    ress2 = k1resn[42,:]

    resv3 = k1resr[1,:]
    ress3 = k1resr[2,:]

    resv4 = k1resr[41,:]
    ress4 = k1resr[42,:]

    errorbar(d, [resv1, resv2, resv3, resv4], [ress1, ress2, ress3, ress4],
             xlabel="$r/a$",
             ylabel="$\\langle W\\rangle$",
             markers=allmarkers,
             lineWidth=0.5,
             legends=["$\\langle W_{s.o.}\\rangle$ at $\\beta=5.65$","$\\langle W_{s.o.}\\rangle$ at $\\beta=6.15$","$\\langle W_{o.o.}\\rangle$ at $\\beta=5.65$","$\\langle W_{o.o.}\\rangle$ at $\\beta=6.15$"],
             savefile=f"../data/distance/BSQK{k+1}.pdf")


fk1resn = np.loadtxt(f'{folderstart}\\WilsonLoopRes\\BS01DistanceFoldN.csv', delimiter=',')
fk1resr = np.loadtxt(f'{folderstart}\\WilsonLoopRes\\BS01DistanceFoldR.csv', delimiter=',')
for k in range(4):
    d = fk1resn[0,:]
    resv1 = fk1resn[1 + 4 * k,:]
    ress1 = fk1resn[2 + 4 * k,:]

    resv2 = fk1resn[240 + 1 + 4 * k,:]
    ress2 = fk1resn[240 + 2 + 4 * k,:]

    resv3 = fk1resr[1 + 4 * k,:]
    ress3 = fk1resr[2 + 4 * k,:]

    resv4 = fk1resr[240 + 1 + 4 * k,:]
    ress4 = fk1resr[240 + 2 + 4 * k,:]

    errorbar(d, [resv1, resv2, resv3, resv4], [ress1, ress2, ress3, ress4],
             xlabel="$r/a$",
             ylabel="$\\langle W\\rangle$",
             markers=allmarkers,
             lineWidth=0.5,
             legends=["$\\langle W_{s.o.}\\rangle$ at $\\beta=5.3$","$\\langle W_{s.o.}\\rangle$ at $\\beta=5.8$",
                      "$\\langle W_{o.o.}\\rangle$ at $\\beta=5.3$","$\\langle W_{o.o.}\\rangle$ at $\\beta=5.8$"],
             savefile=f"../data/distance/BS01FoldK{k+1}.pdf")

fk1resn = np.loadtxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceFoldN.csv', delimiter=',')
fk1resr = np.loadtxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceFoldR.csv', delimiter=',')
for k in range(4):
    d = fk1resn[0, :]
    resv1 = fk1resn[1 + 4 * k, :]
    ress1 = fk1resn[2 + 4 * k, :]

    resv2 = fk1resn[240 + 1 + 4 * k, :]
    ress2 = fk1resn[240 + 2 + 4 * k, :]

    resv3 = fk1resr[1 + 4 * k, :]
    ress3 = fk1resr[2 + 4 * k, :]

    resv4 = fk1resr[240 + 1 + 4 * k, :]
    ress4 = fk1resr[240 + 2 + 4 * k, :]

    errorbar(d, [resv1, resv2, resv3, resv4], [ress1, ress2, ress3, ress4],
             xlabel="$r/a$",
             ylabel="$\\langle W\\rangle$",
             markers=allmarkers,
             lineWidth=0.5,
             legends=["$\\langle W_{s.o.}\\rangle$ at $\\beta=5.65$", "$\\langle W_{s.o.}\\rangle$ at $\\beta=6.15$",
                      "$\\langle W_{o.o.}\\rangle$ at $\\beta=5.65$", "$\\langle W_{o.o.}\\rangle$ at $\\beta=6.15$"],
             savefile=f"../data/distance/BSQFoldK{k + 1}.pdf")
