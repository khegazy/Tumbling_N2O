import glob
import numpy as np
#import scipy.special as sp
import matplotlib.pyplot as plt


def get_bases(folderName, Nbases):
  bases = {}
  files = "{}/*.dat".format(folderName)
  for fl in glob.glob(files):
    L = int(fl[fl.find("_L-")+3:fl.find("_M-")])
    with open(fl, "rb") as file:
      bases[L] = np.fromfile(file, np.double)
      if L != 0:
        bases[L] -= np.mean(bases[L])
      bases[L] /= np.sqrt(np.sum(bases[L]**2))

  basisList = []
  for l in 2*np.arange(Nbases):
    basisList.append(np.expand_dims(bases[l], axis=0))

  return np.concatenate(basisList, axis=0)


def get_data(folder, times):
  
  data = []
  for tm in times:
    tm = np.round(tm*1e3)/1e3
    fName = folder + "molDiff_time-{}_bins[1024,1024].dat".format(tm)
    molDiff = np.fromfile(fName, np.double)
    fName = folder + "atmDiff_time-{}_bins[1024,1024].dat".format(tm)
    atmDiff = np.fromfile(fName, np.double)
    data.append(np.expand_dims(molDiff/atmDiff, axis=0))
  return np.concatenate(data, axis=0)


  
dTs = np.arange(40, dtype=float)*0.1
times = dTs + 316.114

#####  Import Bases/Data  #####
bases = get_bases("/reg/d/psdm/amo/amoi0314/scratch/UED_tumbling/YlmExpVals/UED", 2)

data = get_data(
    "/reg/d/psdm/amo/amoi0314/scratch/diffPatterns/", times)

X,Y   = np.meshgrid(np.arange(-512, 512), np.arange(-512, 512))
print("XY",X,Y)
ratio = X/Y
ratio[512,512] = 0
angs  = np.arctan(ratio)
# Change for X/Y vs Y/X
angs = np.abs(angs)
angs[:512,:] = np.pi - angs[:512,:]
print("nans", np.where(np.isnan(angs)))
R     = np.sqrt(X**2 + Y**2)
print("angs", angs)
plt.pcolormesh(X,Y,angs)
plt.colorbar()
plt.savefig("testAngs.png")
plt.close()
R = np.reshape(R, (-1))
angs = np.reshape(angs, (-1))


cosT = np.cos(angs)
plt.pcolormesh(X,Y,np.reshape(cosT, (1024,1024)))
plt.colorbar()
plt.savefig("testCosT.png")
plt.close()

inds = np.where((R>200) & (R<205))[0]
plt.plot(angs[inds], cosT[inds],'.')

plt.savefig("testCosT_LO.png")
plt.close()

L0 = np.ones_like(R)
#L2P = sp.legendre(2)
#print("poly coeff", L2P)
#print("poly coeff", L2P[2], L2P[1], L2P[0])
L2 = 1.5*cosT**2 - 0.5

plt.pcolormesh(X,Y,np.reshape(L2, (1024,1024)))
plt.colorbar()
plt.savefig("testL2.png")
plt.close()

data_L0 = data*np.expand_dims(L0, axis=0)
data_L2 = data*np.expand_dims(L2, axis=0)
L0_I    = np.zeros((data.shape[0], 256))
L2_I    = np.zeros((data.shape[0], 256))

for i,r in enumerate(2*np.arange(256)):
  inds = np.where((R>=r) & (R<r+2))[0]
  #print("testing jac", r, 0.5*np.sum((L2[inds]**2)/R[inds], axis=-1)/np.pi)
  L0_I[:,i] = 0.5*np.sum(data_L0[:,inds]/R[inds], axis=-1)/np.pi
  L2_I[:,i] = 0.5*np.sum(data_L2[:,inds]/R[inds], axis=-1)/np.pi

L0_I -= np.mean(L0_I, axis=0)
L2_I -= np.mean(L2_I, axis=0)

T,R = np.meshgrid(times, np.arange(256)*2)
plt.pcolormesh(T,R,L0_I.T, vmin=-0.5, vmax=0.5, cmap='seismic')
plt.xlim([times[0], times[-1]])
plt.ylim([0, 512])
plt.colorbar()
plt.savefig("projection_L0.png")
plt.close()

plt.pcolormesh(T,R,L2_I.T, vmin=-0.1, vmax=0.1, cmap='seismic')
plt.xlim([times[0], times[-1]])
plt.ylim([0, 512])
plt.colorbar()
plt.savefig("projection_L2.png")
plt.close()


for i in range(50):
  print("plotting", i)
  tdlo_2 = L2_I[:,4*i]
  rat = tdlo_2[27]/bases[1][27]
  plt.plot(tdlo_2, '-k')
  plt.plot(bases[1]*rat, '-b')
  plt.savefig("compareBasis2_lo-{}.png".format(4*i))
  plt.close()
  
for i in range(50):
  print("plotting0", i)
  tdlo_2 = L2_I[:,4*i]
  rat = tdlo_2[27]/bases[1][27]
  plt.plot(tdlo_2, '-k')
  plt.plot(bases[1]*rat, '-b')
plt.savefig("compareBasis2_lo.png")
