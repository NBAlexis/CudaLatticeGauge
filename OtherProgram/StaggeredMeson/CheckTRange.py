import matplotlib.pyplot as plt

from UsefulFunctions import fitcorrelation

nt = 48
pathstr = "./data/c32p31/correlationp2p"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

for channel in range(20):
    tlst = []
    mplst = []
    mpelst = []
    mmlst = []
    mmelst = []
    chidoflst = []
    for t_start in range(15):
        tlst.append(t_start)
        mp, mpe, mm, mme, chidof = fitcorrelation(channel, pathstr, nt, t_start)
        mplst.append(mp)
        mpelst.append(mpe)
        mmlst.append(mm)
        mmelst.append(mme)
        chidoflst.append(chidof)

    ax1.errorbar(tlst, mplst, yerr=mpelst, fmt='-o')
    ax1.errorbar(tlst, mmlst, yerr=mmelst, fmt=':o')
    ax2.plot(tlst, chidoflst)
ax1.set_xlabel("$t_{start}$")
ax1.set_ylabel("$am$")
ax2.set_xlabel("$t_{start}$")
ax2.set_ylabel("$\chi/d.o.f.$")
plt.show()



