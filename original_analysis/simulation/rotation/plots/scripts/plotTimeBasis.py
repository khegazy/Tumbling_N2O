import numpy as np
import glob
import sys, os
import matplotlib.pyplot as plt

baseDir = "/reg/d/psdm/amo/amoi0314/scratch/UED_tumbling/YlmExpVals/"


for fld in ["UED"]:#, "LCLS", "exp1", "exp2"]:
  files = glob.glob(os.path.join(baseDir, fld, "expValYlm_*"))
  Nfiles = len(files)
  bases = {}
  yTicks = []
  fig, ax = plt.subplots()
  for fln in files:
    print(fln)
    L = int(fln[fln.find("L-")+2:fln.find("_M-")])
    print(L)
    A = np.fromfile(fln, dtype=np.double)
    A /= np.sqrt(np.sum(A**2))
    if fld == "UED":
      A *= 2
    else:
      A *= 5
    """
    A -= np.mean(A)
    if L:
      A /= np.amax(np.abs(A))
      print(np.sum(A**2))
      if np.sum(A**2):
        for i in range(5):
          A /= np.sum(A**2)
      print(np.sum(A**2))
    """
    bases[L] = A
    ax.plot(A+L, '-k')
    yTicks.append(L)
    Nbins = len(A) - 1
  ax.set_ylim([-1, 2*Nfiles - 0.5])
  ax.set_ylabel("L")
  ax.set_xlabel("Time [arb]")
  ax.set_yticks(yTicks)
  ax.set_xlim([0,Nbins])
  fig.tight_layout()
  fig.savefig("../timeBasis_{}.png".format(fld))


  """
  for L,b in bases.iteritems():
    for Ll,bb in bases.iteritems():
      if Ll == L:
        continue
      print(L,Ll,np.sum(b*bb))
  """
