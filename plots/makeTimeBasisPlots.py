import sys
sys.argv.append("-b")
import ROOT 
import numpy as np
import pickle as pl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = u'/reg/neh/home5/khegazy/analysis/tumblingN2O/movie/env/bin/ffmpeg'

#import root_numpy

#from root_numpy import hist2array


sys.path.insert(0, "/reg/neh/home/khegazy/baseTools/tools/")
from pyPlotFunctions import PLOTclass 
#from plotGstyle import genStyle, gen2dStyle


plc = PLOTclass()
ROOT.gROOT.SetStyle("genStyle")
xpix = 1500
ypix = 750
ByPad = 0.05
BxPad = 0.03
BxSep = 0.01
ypixRat = (1-2*ByPad)/4
xpixRat = 0.75
canComp = ROOT.TCanvas("canComp", "canComp", xpix, ypix);

canComp.SetRightMargin(0.01)
canComp.SetLeftMargin(0.085)
canComp.SetTopMargin(0.11)
canComp.SetBottomMargin(0.11)

NfftPows = 15
NfftDers = 2
NfftBasis = NfftPows*NfftDers
ufftBasis = np.fromfile("../rotorBasis/uedfftRotorbasisFine.dat", dtype = np.double)
ufftBasis = np.reshape(ufftBasis, [NfftBasis, -1])
lfftBasis = np.fromfile("../rotorBasis/lclsfftRotorbasis.dat", dtype = np.double)
lfftBasis = np.reshape(lfftBasis, [NfftBasis, -1])
#for i in range(NfftBasis):
#  print(ufftBasis[i])
#print(ufftBasis[1])
#plt.plot(ufftBasis[0,:])
#plt.show()
#plt.plot(ufftBasis[1,:])
#plt.show()

uCosEV = np.fromfile("../simulation/rotation/output/expVals/expValCos2_315.914000-320.214000-0.010000.dat", dtype = np.double)
lCosEV = np.fromfile("/reg/d/psdm/amo/amoi0314/scratch/cosExpVals/expValCos2_317.000000-320.500000-0.004230.dat", dtype=np.double)
uIndMax = np.argmax(uCosEV)
uIndMin = np.argmin(uCosEV)
uIndRef = uIndMin + np.argmin(np.abs(uCosEV[uIndMin:uIndMax]-0.333333))
lIndMax = np.argmax(lCosEV)
lIndMin = np.argmin(lCosEV[:lIndMax])
lIndRef = lIndMin + np.argmin(np.abs(lCosEV[lIndMin:lIndMax]-0.333333))


"""
tmp = np.fromfile("../movie/data/expValCos2_355.726000-360.026000.dat", dtype = np.double)
uCosEV = np.zeros((6,tmp.shape[0]), dtype = np.double)
uCosEV[0,:] = tmp
for i in range(1, uCosEV.shape[0]):
  uCosEV[i,:] = np.fromfile("../movie/data/expValCos" + str((i+1)*2) 
      + "_355.726000-360.026000.dat", dtype = np.double)
lCosEV = np.fromfile("/reg/d/psdm/amo/amoi0314/scratch/cosExpVals/expValCos2_357.045000-360.025000.dat", dtype = np.double)
print("lcosSize",lCosEV.shape)
"""

Nlg = 4
NuRad = 70
itr = 0
ufftCoefRad = np.zeros((Nlg, NfftBasis, NuRad), dtype=np.double) #[Nlg x Ncos x NuRad]
for ilg in range(Nlg):
  itr = 0
  for ip in range(NfftPows):
    for idr in range(NfftDers):
      ufftCoefRad[ilg,itr,:]  = np.fromfile("../movie/data/UEDfft_" + str(ip) 
          + "-" + str(idr) + "_CoefLeg" + str(ilg*2) 
          + "rad" + str(NuRad) + ".dat", dtype = np.double)
      #print(ilg, ip, idr)
      #if ilg == 2:
      #  print(ufftCoefRad[ilg,itr,:])
      itr += 1

NlRad = 30
#NlRad = 60
# Sum bin 10 - 14
lfftCoefRad = np.zeros((Nlg, NfftBasis, NlRad), dtype=np.double)
for ilg in range(Nlg):
  #if (ilg>2):
  #  continue
  lInp = np.loadtxt("../movie/data/projections-amoi0314-r0172_e5_g1_l" + str(2*ilg) 
     + "_tspan1.083.dat", dtype = np.double)
  #lInp1 = np.loadtxt("../movie/data/projections-amoi0314-r0172_e9_g1_l" + str(2*ilg) 
  #   + "_tspan1.083.dat", dtype = np.double)
  #lInp2 = np.loadtxt("../movie/data/projections-amoi0314-r0172_e10_g1_l" + str(2*ilg) 
  #   + "_tspan1.083.dat", dtype = np.double)
  for i in range(NfftBasis):
    lfftCoefRad[ilg,i,:] = lInp[:,i+1]
  #for i in range(NfftBasis):
  #  lfftCoefRad[ilg,i,:] = lInp1[:,i+1] + lInp2[:,i+1]
    



####################
###  Make Plot  ###
####################

maxQ = 14.1207
Npix = 257
circRad = 0.05
dX = np.tile(np.arange(Npix) - (Npix - 1)/2, (Npix,1))
dR = np.sqrt(dX**2 + np.transpose(dX)**2)
dR[dR==0] = 1
cosTheta = dX/dR

#print(dR[Npix/2+1,:])

# radial indices
ufftCoef = np.zeros((Nlg, NfftBasis, Npix, Npix))
rInds = np.copy(dR)
rInds[rInds > Npix/2] = 0
rInds *= 1.0*NuRad/(1.0*Npix/2.)
rIndsFlat_double = np.reshape(rInds, (-1, 1))
rIndsFlat = rIndsFlat_double.astype(np.int8)
for ilg in range(Nlg):
  for i in range(NfftBasis):
    ufftCoef[ilg,i,:,:] = np.reshape(ufftCoefRad[ilg,i,rIndsFlat], (Npix, Npix))


lfftCoef = np.zeros((Nlg, NfftBasis, Npix, Npix))
rInds = np.copy(dR)
rInds[rInds > Npix/2] = 0
rInds -= circRad*Npix/2.
rInds[rInds < 0] = 0
rInds *= 1.0*NlRad/((1 - circRad)*Npix/2.)
rIndsFlat_double = np.reshape(rInds, (-1, 1))
rIndsFlat = rIndsFlat_double.astype(np.int8)
for ilg in range(Nlg):
  for i in range(NfftBasis):
    lfftCoef[ilg,i,:,:] = np.reshape(lfftCoefRad[ilg,i,rIndsFlat], (Npix, Npix))

#plt.imshow(rInds)
#plt.show()


legImg = np.zeros((Nlg, Npix, Npix), dtype = np.double)
for ilg in range(Nlg):
  c = np.zeros(2*ilg + 1)
  c[-1] = 1
  legImg[ilg,:,:] = np.polynomial.legendre.legval(cosTheta, c)

"""
img = np.zeros((Npix, Npix), dtype = np.double)
for ilg in range(1, Nlg):
  for ibs in range(NfftBasis):
    img[:,:] += ufftCoef[ilg,ibs,:,:]*ufftBasis[ibs,25]*legImg[ilg,:,:]
"""


##################
###  Make Plot  ##
##################

qmin = 0; qmax = 1;
uminHist = ROOT.TH2F("uminHist", "uminHist", Npix/2+1, qmin, qmax, Npix, qmin, qmax)
lminHist = ROOT.TH2F("lminHist", "lminHist", Npix/2+1, qmin, qmax, Npix, qmin, qmax)
umaxHist = ROOT.TH2F("umaxHist", "umaxHist", Npix/2+1, qmin, qmax, Npix, qmin, qmax)
lmaxHist = ROOT.TH2F("lmaxHist", "lmaxHist", Npix/2+1, qmin, qmax, Npix, qmin, qmax)
urefHist = ROOT.TH2F("urefHist", "urefHist", Npix/2+1, qmin, qmax, Npix, qmin, qmax)
lrefHist = ROOT.TH2F("lrefHist", "lrefHist", Npix/2+1, qmin, qmax, Npix, qmin, qmax)

uminHist.SetContour(80)
lminHist.SetContour(80)
umaxHist.SetContour(80)
lmaxHist.SetContour(80)
urefHist.SetContour(80)
lrefHist.SetContour(80)

# Reference Image
uRefImg = np.zeros((Npix, Npix), dtype = np.double)
lRefImg = np.zeros((Npix, Npix), dtype = np.double)
#for ilg in range(1,Nlg):
for ilg in range(Nlg):
  for ibs in range(NfftBasis):
    lRefImg[:,:] += lfftCoef[ilg,ibs]*ufftBasis[ibs,uIndRef]*legImg[ilg,:,:]
    if ilg is not 0:
      uRefImg[:,:] += ufftCoef[ilg,ibs,:,:]*ufftBasis[ibs,uIndRef]*legImg[ilg,:,:]


# Minimum Image
uMinImg = np.zeros((Npix, Npix), dtype = np.double)
lMinImg = np.zeros((Npix, Npix), dtype = np.double)
#for ilg in range(1,Nlg):
for ilg in range(Nlg):
  for ibs in range(NfftBasis):
    lMinImg[:,:] += lfftCoef[ilg,ibs]*ufftBasis[ibs,uIndMin]*legImg[ilg,:,:]
    if ilg is not 0:
      uMinImg[:,:] += ufftCoef[ilg,ibs,:,:]*ufftBasis[ibs,uIndMin]*legImg[ilg,:,:]

# Maximum Image
uMaxImg = np.zeros((Npix, Npix), dtype = np.double)
lMaxImg = np.zeros((Npix, Npix), dtype = np.double)
#for ilg in range(1,Nlg):
for ilg in range(Nlg):
  for ibs in range(NfftBasis):
    lMaxImg[:,:] += lfftCoef[ilg,ibs]*ufftBasis[ibs,uIndMax]*legImg[ilg,:,:]
    if ilg is not 0:
      uMaxImg[:,:] += ufftCoef[ilg,ibs,:,:]*ufftBasis[ibs,uIndMax]*legImg[ilg,:,:]

uNorm = 1.0/max(np.absolute([np.amax(uRefImg), np.amax(uMinImg), np.amax(uMaxImg), np.amin(uRefImg), np.amin(uMinImg), np.amin(uMaxImg)]))

lNorm = 1.0/max(np.absolute([np.amax(lRefImg), np.amax(lMinImg), np.amax(lMaxImg), np.amin(lRefImg), np.amin(lMinImg), np.amin(lMaxImg)]))

uMinImg *= uNorm
uRefImg *= uNorm
uMaxImg *= uNorm
lMinImg *= lNorm
lRefImg *= lNorm
lMaxImg *= lNorm

for ir in range(Npix):
  for ic in range(Npix/2 + 1):
    lrefHist.SetBinContent(ic+1, ir+1, lRefImg[ir,ic])
    lminHist.SetBinContent(ic+1, ir+1, lMinImg[ir,ic])
    lmaxHist.SetBinContent(ic+1, ir+1, lMaxImg[ir,ic])
    urefHist.SetBinContent(ic+1, ir+1, uRefImg[ir,Npix/2 + ic])
    uminHist.SetBinContent(ic+1, ir+1, uMinImg[ir,Npix/2 + ic])
    umaxHist.SetBinContent(ic+1, ir+1, uMaxImg[ir,Npix/2 + ic])


#norm = 1.0/max(np.absolute([refHist.GetMaximum(), minHist.GetMaximum(), maxHist.GetMaximum(), refHist.GetMinimum(), minHist.GetMinimum(), maxHist.GetMinimum()]))

for ir in range(Npix):
  for ic in range(Npix/2 + 1):
    if np.sqrt((ir-Npix/2)**2 + (ic-Npix/2)**2) > Npix/2:
      umaxHist.SetBinContent(Npix/2 - ic + 2, ir+1, -1e20)
      uminHist.SetBinContent(Npix/2 - ic + 2, ir+1, -1e20)
      urefHist.SetBinContent(Npix/2 - ic + 2, ir+1, -1e20)
      lmaxHist.SetBinContent(ic+1, ir+1, -1e20)
      lminHist.SetBinContent(ic+1, ir+1, -1e20)
      lrefHist.SetBinContent(ic+1, ir+1, -1e20)


#refHist.Scale(norm)
#maxHist.Scale(norm)
#minHist.Scale(norm)
urefHist.SetMaximum(1)
umaxHist.SetMaximum(1)
uminHist.SetMaximum(1)
lrefHist.SetMaximum(1)
lmaxHist.SetMaximum(1)
lminHist.SetMaximum(1)
urefHist.SetMinimum(-1)
umaxHist.SetMinimum(-1)
uminHist.SetMinimum(-1)
lrefHist.SetMinimum(-1)
lmaxHist.SetMinimum(-1)
lminHist.SetMinimum(-1)

uTi = 316 #355.726
uTf = 320.5 #360.026
timeRatio = 1.95
timeWindow = 150
lOrigTi = 317.0
lOrigTf = 320.5
udt = (uTf - uTi)/ufftBasis[0].shape[0]
ldt = (lOrigTf - lOrigTi)/lfftBasis.shape[1]
uTref = udt*uIndRef + uTi
lTref = ldt*lIndRef + lOrigTi
lShift = 2*timeRatio*(uTref - (lOrigTi - (lTref - uTref)))*(lIndRef - lfftBasis.shape[1]/2.)/float(lfftBasis.shape[1])
print("Time of full revival UED/LCLS: ",uTref, lTref)
uCosHist = ROOT.TH1F("UEDcosSq", "UEDcosSq", \
    ufftBasis[0].shape[0], uTi, uTf)
lCosHist = ROOT.TH1F("LCLScosSq", "LCLScosSq", \
    lfftBasis.shape[1], 
    lShift + uTref - timeRatio*(uTref - (lOrigTi - (lTref - uTref))), 
    lShift + uTref + timeRatio*((lOrigTf - (lTref - uTref) - uTref)))
uB1Hist = ROOT.TH1F("UEDb1", "UEDb1", \
    ufftBasis[0].shape[0], uTi, uTf)
lB1Hist = ROOT.TH1F("LCLSb1", "LCLSb1", \
    lfftBasis.shape[1], 
    lShift + uTref - timeRatio*(uTref - (lOrigTi - (lTref - uTref))), 
    lShift + uTref + timeRatio*((lOrigTf - (lTref - uTref) - uTref)))
uB2Hist = ROOT.TH1F("UEDb2", "UEDb2", \
    ufftBasis[0].shape[0], uTi, uTf)
lB2Hist = ROOT.TH1F("LCLSb2", "LCLSb2", \
    lfftBasis.shape[1], 
    lShift + uTref - timeRatio*(uTref - (lOrigTi - (lTref - uTref))), 
    lShift + uTref + timeRatio*((lOrigTf - (lTref - uTref) - uTref)))
uB3Hist = ROOT.TH1F("UEDb3", "UEDb3", \
    ufftBasis[0].shape[0], uTi, uTf)
lB3Hist = ROOT.TH1F("LCLSb3", "LCLSb3", \
    lfftBasis.shape[1], 
    lShift + uTref - timeRatio*(uTref - (lOrigTi - (lTref - uTref))), 
    lShift + uTref + timeRatio*((lOrigTf - (lTref - uTref) - uTref)))
uB4Hist = ROOT.TH1F("UEDb4", "UEDb4", \
    ufftBasis[0].shape[0], uTi, uTf)
lB4Hist = ROOT.TH1F("LCLSb4", "LCLSb4", \
    lfftBasis.shape[1], 
    lShift + uTref - timeRatio*(uTref - (lOrigTi - (lTref - uTref))), 
    lShift + uTref + timeRatio*((lOrigTf - (lTref - uTref) - uTref)))
uB5Hist = ROOT.TH1F("UEDb5", "UEDb5", \
    ufftBasis[0].shape[0], uTi, uTf)
lB5Hist = ROOT.TH1F("LCLSb5", "LCLSb5", \
    lfftBasis.shape[1], 
    lShift + uTref - timeRatio*(uTref - (lOrigTi - (lTref - uTref))), 
    lShift + uTref + timeRatio*((lOrigTf - (lTref - uTref) - uTref)))
uB6Hist = ROOT.TH1F("UEDb6", "UEDb6", \
    ufftBasis[0].shape[0], uTi, uTf)
lB6Hist = ROOT.TH1F("LCLSb6", "LCLSb6", \
    lfftBasis.shape[1], 
    lShift + uTref - timeRatio*(uTref - (lOrigTi - (lTref - uTref))), 
    lShift + uTref + timeRatio*((lOrigTf - (lTref - uTref) - uTref)))
uB7Hist = ROOT.TH1F("UEDb7", "UEDb7", \
    ufftBasis[0].shape[0], uTi, uTf)
lB7Hist = ROOT.TH1F("LCLSb7", "LCLSb7", \
    lfftBasis.shape[1], 
    lShift + uTref - timeRatio*(uTref - (lOrigTi - (lTref - uTref))), 
    lShift + uTref + timeRatio*((lOrigTf - (lTref - uTref) - uTref)))








for i in range(ufftBasis[0].shape[0]):
  uCosHist.SetBinContent(i+1, ufftBasis[0,i])
  uB1Hist.SetBinContent(i+1, ufftBasis[2,i])
  uB2Hist.SetBinContent(i+1, ufftBasis[4,i])
  uB3Hist.SetBinContent(i+1, ufftBasis[6,i])
  uB4Hist.SetBinContent(i+1, ufftBasis[8,i])
  uB5Hist.SetBinContent(i+1, ufftBasis[10,i])
  uB6Hist.SetBinContent(i+1, ufftBasis[12,i])
  uB7Hist.SetBinContent(i+1, ufftBasis[14,i])

for i in range(lfftBasis.shape[1]):
  lCosHist.SetBinContent(i+1, lfftBasis[0,i])
  lB1Hist.SetBinContent(i+1, lfftBasis[2,i])
  lB2Hist.SetBinContent(i+1, lfftBasis[4,i])
  lB3Hist.SetBinContent(i+1, lfftBasis[6,i])
  lB4Hist.SetBinContent(i+1, lfftBasis[8,i])
  lB5Hist.SetBinContent(i+1, lfftBasis[10,i])
  lB6Hist.SetBinContent(i+1, lfftBasis[12,i])
  lB7Hist.SetBinContent(i+1, lfftBasis[14,i])

print("limits", timeWindow*udt, uTref - timeWindow*udt, uTref + timeWindow*udt)
uCosHist.GetXaxis().SetRangeUser(uTref - timeWindow*udt, uTref + timeWindow*udt)
uB1Hist.GetXaxis().SetRangeUser(uTref - timeWindow*udt, uTref + timeWindow*udt)
uB2Hist.GetXaxis().SetRangeUser(uTref - timeWindow*udt, uTref + timeWindow*udt)
uB3Hist.GetXaxis().SetRangeUser(uTref - timeWindow*udt, uTref + timeWindow*udt)
uB4Hist.GetXaxis().SetRangeUser(uTref - timeWindow*udt, uTref + timeWindow*udt)
uB5Hist.GetXaxis().SetRangeUser(uTref - timeWindow*udt, uTref + timeWindow*udt)
uB6Hist.GetXaxis().SetRangeUser(uTref - timeWindow*udt, uTref + timeWindow*udt)
uB7Hist.GetXaxis().SetRangeUser(uTref - timeWindow*udt, uTref + timeWindow*udt)


UEDfitFile = ROOT.TFile.Open("../simulation/fitParameters/UEDchiSqVtimeFit.root","read")
UEDfit = UEDfitFile.Get("fit")
UEDfit.GetXaxis().SetRangeUser(-83, 85);
LCLSfitFile = ROOT.TFile.Open("../simulation/fitParameters/LCLSchiSqVtimeFit.root","read")
LCLSfit = LCLSfitFile.Get("fit")
LCLSfit.GetXaxis().SetRangeUser(-30, 30);
uColor = 4
lColor = 2
labelSize = 0.04
UEDfit.SetMarkerSize(0.6)
LCLSfit.SetMarkerSize(0.6)
UEDfit.SetMarkerColor(uColor)
uCosHist.SetLineColor(uColor)
uB1Hist.SetLineColor(uColor)
uB2Hist.SetLineColor(uColor)
uB3Hist.SetLineColor(uColor)
uB4Hist.SetLineColor(uColor)
uB5Hist.SetLineColor(uColor)
uB6Hist.SetLineColor(uColor)
uB7Hist.SetLineColor(uColor)
LCLSfit.SetMarkerColor(lColor)
lCosHist.SetLineColor(lColor)
lB1Hist.SetLineColor(lColor)
lB2Hist.SetLineColor(lColor)
lB3Hist.SetLineColor(lColor)
lB4Hist.SetLineColor(lColor)
lB5Hist.SetLineColor(lColor)
lB6Hist.SetLineColor(lColor)
lB7Hist.SetLineColor(lColor)
uCosHist.GetXaxis().SetTitleColor(uColor)
uB1Hist.GetXaxis().SetTitleColor(uColor)
uB2Hist.GetXaxis().SetTitleColor(uColor)
uB3Hist.GetXaxis().SetTitleColor(uColor)
uB4Hist.GetXaxis().SetTitleColor(uColor)
uB5Hist.GetXaxis().SetTitleColor(uColor)
uB6Hist.GetXaxis().SetTitleColor(uColor)
uB7Hist.GetXaxis().SetTitleColor(uColor)
uCosHist.GetXaxis().SetLabelColor(uColor)
#uB1Hist.GetYaxis().SetTitle("#LT Cos^{2}(#theta) #GT")
#uB2Hist.GetYaxis().SetTitle("#LT Cos^{2}(#theta) #GT")
uCosHist.GetYaxis().SetTitle("[Arb]")##LT Cos^{2}(#theta) #GT")
uB1Hist.GetYaxis().SetTitle("[Arb]")#LT Cos^{2}(#theta) #GT")
uB2Hist.GetYaxis().SetTitle("[Arb]")#LT Cos^{2}(#theta) #GT")
uB3Hist.GetYaxis().SetTitle("[Arb]")#LT Cos^{2}(#theta) #GT")
uB4Hist.GetYaxis().SetTitle("[Arb]")#LT Cos^{2}(#theta) #GT")
uB5Hist.GetYaxis().SetTitle("[Arb]")#LT Cos^{2}(#theta) #GT")
uB6Hist.GetYaxis().SetTitle("[Arb]")#LT Cos^{2}(#theta) #GT")
uB7Hist.GetYaxis().SetTitle("[Arb]")#LT Cos^{2}(#theta) #GT")
uB7Hist.GetXaxis().SetTitle("Time [ps]")
uB6Hist.GetXaxis().SetTitle("Time [ps]")
UEDfit.GetXaxis().CenterTitle(1)
LCLSfit.GetXaxis().CenterTitle(1)
uB6Hist.GetXaxis().CenterTitle(1)
uB7Hist.GetXaxis().CenterTitle(1)
UEDfit.GetYaxis().CenterTitle(1)
LCLSfit.GetYaxis().CenterTitle(1)
uCosHist.GetYaxis().CenterTitle(1)
uB1Hist.GetYaxis().CenterTitle(1)
uB2Hist.GetYaxis().CenterTitle(1)
uB3Hist.GetYaxis().CenterTitle(1)
uB4Hist.GetYaxis().CenterTitle(1)
uB5Hist.GetYaxis().CenterTitle(1)
uB6Hist.GetYaxis().CenterTitle(1)
uB7Hist.GetYaxis().CenterTitle(1)
uCosHist.SetMaximum(1.1*lCosHist.GetMaximum())
uB1Hist.SetMinimum(1.1*lB1Hist.GetMinimum())
uB2Hist.SetMaximum(1.1*lB2Hist.GetMaximum())
uB3Hist.SetMaximum(1.1*lB3Hist.GetMaximum())
uB4Hist.SetMaximum(1.1*lB4Hist.GetMaximum())
uB5Hist.SetMaximum(1.1*lB5Hist.GetMaximum())
uB6Hist.SetMaximum(1.1*lB6Hist.GetMaximum())
uB7Hist.SetMaximum(1.1*lB7Hist.GetMaximum())
uCosHist.GetYaxis().SetNdivisions(503)
uB1Hist.GetYaxis().SetNdivisions(503)
uB2Hist.GetYaxis().SetNdivisions(503)
uB3Hist.GetYaxis().SetNdivisions(503)
uB4Hist.GetYaxis().SetNdivisions(503)
uB5Hist.GetYaxis().SetNdivisions(503)
uB6Hist.GetYaxis().SetNdivisions(503)
uB7Hist.GetYaxis().SetNdivisions(503)

umaxHist.GetXaxis().SetLabelSize(0)
umaxHist.GetYaxis().SetLabelSize(0)
umaxHist.GetXaxis().SetTickSize(0)
umaxHist.GetYaxis().SetTickSize(0)
umaxHist.GetYaxis().SetTitleColor(0)
lmaxHist.GetXaxis().SetLabelSize(0)
lmaxHist.GetYaxis().SetLabelSize(0)
lmaxHist.GetXaxis().SetTickSize(0)
lmaxHist.GetYaxis().SetTickSize(0)
lmaxHist.GetYaxis().SetTitleColor(0)

uminHist.GetXaxis().SetLabelSize(0)
uminHist.GetXaxis().SetLabelOffset(999)
uminHist.GetYaxis().SetLabelSize(0)
uminHist.GetXaxis().SetTickSize(0)
uminHist.GetYaxis().SetTickSize(0)
lminHist.GetXaxis().SetLabelSize(0)
lminHist.GetXaxis().SetLabelOffset(999)
lminHist.GetYaxis().SetLabelSize(0)
lminHist.GetXaxis().SetTickSize(0)
lminHist.GetYaxis().SetTickSize(0)

"""
uCosHist.Draw()
canComp.Print("utest.png")
lCosHist.Draw()
canComp.Print("ltest.png")
"""
canComp.cd()
padB0 = ROOT.TPad("padB0", "padB0", 0, 1 - ByPad, 1, 1)
padB0.Draw()
padB0.cd()


canComp.cd()
padB5 = ROOT.TPad("padB5", "padB5", 0, 0, 1, ByPad)


#plt.plot(l1dSimAln)
#plt.savefig("l1dSimAln")


txt = ROOT.TText()
tex = ROOT.TLatex()
line = ROOT.TLine()
lb = ROOT.TLine()
rb = ROOT.TLine()
tb = ROOT.TLine()
bb = ROOT.TLine()
lb.SetLineColor(0)
rb.SetLineColor(0)
tb.SetLineColor(0)
bb.SetLineColor(0)
lcirc = ROOT.TEllipse(1.0, 0.5, 2*circRad, circRad, 90, 270)
lcirc.SetFillColor(1)
ucirc = ROOT.TEllipse(0.0, 0.5, 2*circRad, circRad)
ucirc.SetFillColor(1)
uAxis = ROOT.TGaxis(2*circRad, 0.5, 1.0, 0.5, maxQ*circRad, maxQ, 507, "-+")
uAxis.SetLabelSize(labelSize*1.5)
uAxis.SetLabelOffset(-0.008)
uAxis.SetTitleSize(labelSize*1.5)
uAxis.SetTitle("Q [#AA^{-1}]")
lAxis = ROOT.TGaxis(0, 0.5, 1 - 2*circRad, 0.5, 0, 30, 507, "-+")
lAxis.SetLabelSize(0)
laxisLabels = ["530", "525", "520", "515", "510", "505"]
laxisPosItr = (1.0-circRad*2)/(len(laxisLabels))
txt.SetTextSize(labelSize*1.2)
txt.SetTextFont(62)
txt.SetTextColor(1)
txtLabel = ROOT.TText()
txtLabel.SetTextFont(62)
txtLabel.SetTextSize(0.025)


pSize = 0.155
pSep = 0.004
p1x = 0.34
p1y = 0.46 
canComp.cd()
padFU = ROOT.TPad("padFU", "padFU", xpixRat, 0.5 + ByPad/2., 1, 1 - ByPad)
padFU.SetRightMargin(0.05)
padFU.SetLeftMargin(0.15)
padFU.SetTopMargin(0)
padFU.SetBottomMargin(0)
padFU.SetFrameFillColor(0)
padFU.SetFrameLineColor(0)
padFU.SetBorderSize(0)
padFU.Draw()
padFU.cd()
UEDfit.GetXaxis().SetLabelFont(62)
UEDfit.GetXaxis().SetTitleFont(62)
UEDfit.GetYaxis().SetLabelFont(62)
UEDfit.GetYaxis().SetTitleFont(62)
UEDfit.GetXaxis().SetLabelSize(labelSize*1.35)
UEDfit.GetYaxis().SetLabelSize(labelSize*1.35)
UEDfit.GetXaxis().SetTitleSize(labelSize*1.45)
UEDfit.GetYaxis().SetTitleSize(labelSize*1.45)
UEDfit.GetXaxis().SetTitleOffset(0.85)
UEDfit.GetYaxis().SetTitleOffset(0.98)
UEDfit.Draw()
tex.SetTextFont(62)
tex.SetTextSize(0.07)
tex.DrawLatexNDC(0.39,0.85,"UED t_{ref} Fit")


p2x = 0.585
p2y = 0.19
canComp.cd()
padFL = ROOT.TPad("padFL", "padFL", xpixRat, ByPad, 1.0, 0.5 - ByPad/2.)
padFL.SetRightMargin(0.05)
padFL.SetLeftMargin(0.15)
padFL.SetTopMargin(0)
padFL.SetBottomMargin(0)
padFL.SetFrameFillColor(0)
padFL.SetBorderSize(0)
padFL.Draw()
padFL.cd()
LCLSfit.GetXaxis().SetLabelFont(62)
LCLSfit.GetXaxis().SetTitleFont(62)
LCLSfit.GetYaxis().SetLabelFont(62)
LCLSfit.GetYaxis().SetTitleFont(62)
LCLSfit.GetXaxis().SetLabelSize(labelSize*1.35)
LCLSfit.GetYaxis().SetLabelSize(labelSize*1.35)
LCLSfit.GetXaxis().SetTitleSize(labelSize*1.45)
LCLSfit.GetYaxis().SetTitleSize(labelSize*1.45)
LCLSfit.GetXaxis().SetTitleOffset(0.85)
LCLSfit.GetYaxis().SetTitleOffset(0.98)
LCLSfit.Draw()
tex.SetTextFont(62)
tex.SetTextSize(0.07)
tex.DrawLatexNDC(0.39,0.85,"LCLS t_{ref} Fit")

tempC = ROOT.TCanvas("Tcan","Tcan",500,500)
tempC.cd()
tempC.SetLeftMargin(0.18)
tempC.SetBottomMargin(0.18)
LCLSfit.Draw()
tempC.Print("LCLSchiSqFit.png")



canComp.cd()
padX0 = ROOT.TPad("padX0", "padX0", 0, 0, BxPad, 1)
padX0.Draw()
padX0.cd()

Bx1E = BxPad + (xpixRat - BxPad - BxSep)/2
Bx2S = BxPad + (xpixRat - BxPad + BxSep)/2

canComp.cd()
padB1 = ROOT.TPad("padB1", "padB1", BxPad, 0, xpixRat, 1)
padB1.SetRightMargin(0.0)
padB1.SetLeftMargin(0.1)
padB1.SetTopMargin(0.1)
padB1.SetBottomMargin(0.1)
padB1.SetFrameFillColor(0)
padB1.SetFrameLineColor(0)
padB1.SetBorderSize(0)
padB1.Draw()
padB1.cd()
uCosHist.GetXaxis().SetLabelSize(labelSize)
uCosHist.GetXaxis().SetTitleSize(0)
uCosHist.GetXaxis().SetLabelFont(62)
uCosHist.GetXaxis().SetTitleFont(62)
#uCosHist.GetXaxis().SetTitleOffset(0)
uCosHist.GetXaxis().SetLabelColor(uColor)
uCosHist.GetXaxis().SetTitle("Time [ps]")
uCosHist.GetYaxis().SetLabelSize(labelSize)
uCosHist.GetYaxis().SetTitleSize(labelSize)
uCosHist.GetYaxis().SetLabelFont(62)
uCosHist.GetYaxis().SetTitleFont(62)
#uCosHist.GetYaxis().SetTitleOffset(-1.27)

uCosHist.Draw("C")
lCosHist.Draw("SAMEC")

cosLgnd = ROOT.TLegend(0.18, 0.68, 0.32, 0.8)
#cosLgnd = ROOT.TLegend(0.75, 0.7, 0.95, 0.82)
cosLgnd.AddEntry(uCosHist, "UED", "l")
cosLgnd.AddEntry(lCosHist, "LCLS", "l")
cosLgnd.Draw("SAME")

txt.SetTextSize(labelSize)
txt.SetTextAlign(22)
txt.SetTextFont(62)
txt.SetTextColor(lColor)
axLabels = [317.5, 318, 318.5, 319, 319.5, 320]
for i in axLabels:
  txt.DrawText(i, uCosHist.GetMaximum()*1.05, str((i-uTref)/timeRatio + lTref)[:5])
txt.DrawText((uTref-timeWindow*udt + uTref+timeWindow*udt)/2, 
      uCosHist.GetMaximum()*1.27/1.1, "Time [ps]")

txt.SetTextColor(uColor)
txt.DrawText((uTref-timeWindow*udt + uTref+timeWindow*udt)/2, 
      uCosHist.GetMinimum()*1.3, "Time [ps]")
txt.SetTextColor(lColor)

tex.SetTextFont(42)
tex.SetTextSize(0.2)
#tex.DrawLatexNDC(0.05,0.15,"b_{0}")

"""
canComp.cd()
padB2 = ROOT.TPad("padB2", "padB2", Bx2S, 1 - ypixRat - ByPad, xpixRat, 1 - ByPad)
padB2.SetRightMargin(0.0)
padB2.SetLeftMargin(0.0)
padB2.SetTopMargin(0.0)
padB2.SetBottomMargin(0)
padB2.SetFrameFillColor(0)
padB2.SetBorderSize(0)
padB2.Draw()
padB2.cd()
uB1Hist.GetXaxis().SetLabelSize(labelSize*3)
uB1Hist.GetXaxis().SetTitleSize(labelSize*3)
uB1Hist.GetXaxis().SetLabelFont(62)
uB1Hist.GetXaxis().SetTitleFont(62)
uB1Hist.GetXaxis().SetTitleOffset(0.91)
uB1Hist.GetXaxis().SetLabelColor(uColor)
uB1Hist.GetYaxis().SetLabelSize(0)
uB1Hist.GetYaxis().SetTitleSize(0)
uB1Hist.GetYaxis().SetLabelFont(62)
uB1Hist.GetYaxis().SetTitleFont(62)
uB1Hist.GetYaxis().SetTitleOffset(0.27)
uB1Hist.Draw("C")
lB1Hist.Draw("SAMEC")

scl = uCosHist.GetMaximum()/uB1Hist.GetMaximum()
for i in axLabels:
  txt.DrawText(i, uCosHist.GetMaximum()*scl*1.22/1.1, str((i-uTref)/timeRatio + lTref)[:5])
txt.DrawText((uTref-timeWindow*udt + uTref+timeWindow*udt)/2, 
      uCosHist.GetMaximum()*scl*1.44/1.1, "Time [ps]")
tex.SetTextFont(42)
tex.SetTextSize(0.2)
tex.DrawLatexNDC(0.05,0.15,"b_{0}")

canComp.cd()
padB3 = ROOT.TPad("padB3", "padB3", BxPad, 1 - 2*ypixRat - ByPad, Bx1E, 1 - ypixRat - ByPad)
padB3.SetRightMargin(0.0)
padB3.SetLeftMargin(0.0)
padB3.SetTopMargin(0.0)
padB3.SetBottomMargin(0)
padB3.SetFrameFillColor(0)
padB3.SetBorderSize(0)
padB3.Draw()
padB3.cd()
uB2Hist.GetXaxis().SetLabelSize(labelSize*3)
uB2Hist.GetXaxis().SetTitleSize(labelSize*3)
uB2Hist.GetXaxis().SetLabelFont(62)
uB2Hist.GetXaxis().SetTitleFont(62)
uB2Hist.GetXaxis().SetTitleOffset(0.91)
uB2Hist.GetXaxis().SetLabelColor(uColor)
uB2Hist.GetYaxis().SetLabelSize(labelSize*3)
uB2Hist.GetYaxis().SetTitleSize(labelSize*3)
uB2Hist.GetYaxis().SetLabelFont(62)
uB2Hist.GetYaxis().SetTitleFont(62)
uB2Hist.GetYaxis().SetTitleOffset(0.27)
uB2Hist.Draw("C")
lB2Hist.Draw("SAMEC")
#cosLgnd.Draw("SAME")
tex.DrawLatexNDC(0.05,0.15,"b_{1}")
print("here1")

canComp.cd()
padB4 = ROOT.TPad("padB4", "padB4", Bx2S, 1 - 2*ypixRat - ByPad, xpixRat, 1 - ypixRat - ByPad)
padB4.SetRightMargin(0.0)
padB4.SetLeftMargin(0.0)
padB4.SetTopMargin(0.0)
padB4.SetBottomMargin(0.0)
padB4.SetFrameFillColor(0)
padB4.SetBorderSize(0)
padB4.Draw()
padB4.cd()
uB3Hist.GetXaxis().SetLabelSize(labelSize*3)
uB3Hist.GetXaxis().SetTitleSize(labelSize*3)
uB3Hist.GetXaxis().SetLabelFont(62)
uB3Hist.GetXaxis().SetTitleFont(62)
uB3Hist.GetXaxis().SetTitleOffset(0.91)
uB3Hist.GetXaxis().SetLabelColor(uColor)
uB3Hist.GetYaxis().SetLabelSize(0)
uB3Hist.GetYaxis().SetTitleSize(0)
uB3Hist.GetYaxis().SetLabelFont(62)
uB3Hist.GetYaxis().SetTitleFont(62)
uB3Hist.GetYaxis().SetTitleOffset(0.27)
uB3Hist.Draw("C")
lB3Hist.Draw("SAMEC")
tex.DrawLatexNDC(0.05,0.15,"b_{1}^{(1)}")

canComp.cd()
padB5 = ROOT.TPad("padB5", "padB5", BxPad, 1 - 3*ypixRat - ByPad, Bx1E, 1 - 2*ypixRat - ByPad)
padB5.SetRightMargin(0.0)
padB5.SetLeftMargin(0.0)
padB5.SetTopMargin(0.0)
padB5.SetBottomMargin(0.0)
padB5.SetFrameFillColor(0)
padB5.SetBorderSize(0)
padB5.Draw()
padB5.cd()
uB4Hist.GetXaxis().SetLabelSize(labelSize*3)
uB4Hist.GetXaxis().SetTitleSize(labelSize*3)
uB4Hist.GetXaxis().SetLabelFont(62)
uB4Hist.GetXaxis().SetTitleFont(62)
uB4Hist.GetXaxis().SetTitleOffset(0.91)
uB4Hist.GetXaxis().SetLabelColor(uColor)
uB4Hist.GetYaxis().SetLabelSize(labelSize*3)
uB4Hist.GetYaxis().SetTitleSize(labelSize*3)
uB4Hist.GetYaxis().SetLabelFont(62)
uB4Hist.GetYaxis().SetTitleFont(62)
uB4Hist.GetYaxis().SetTitleOffset(0.27)
uB4Hist.Draw("C")
lB4Hist.Draw("SAMEC")
tex.DrawLatexNDC(0.05,0.15,"b_{2}^{(0)}")

canComp.cd()
padB6 = ROOT.TPad("padB6", "padB6", Bx2S, 1 - 3*ypixRat - ByPad, xpixRat, 1 - 2*ypixRat - ByPad)
padB6.SetRightMargin(0.0)
padB6.SetLeftMargin(0.0)
padB6.SetTopMargin(0.0)
padB6.SetBottomMargin(0.0)
padB6.SetFrameFillColor(0)
padB6.SetBorderSize(0)
padB6.Draw()
padB6.cd()
uB5Hist.GetXaxis().SetLabelSize(labelSize*3)
uB5Hist.GetXaxis().SetTitleSize(labelSize*3)
uB5Hist.GetXaxis().SetLabelFont(62)
uB5Hist.GetXaxis().SetTitleFont(62)
uB5Hist.GetXaxis().SetTitleOffset(0.91)
uB5Hist.GetXaxis().SetLabelColor(uColor)
uB5Hist.GetYaxis().SetLabelSize(0)
uB5Hist.GetYaxis().SetTitleSize(0)
uB5Hist.GetYaxis().SetLabelFont(62)
uB5Hist.GetYaxis().SetTitleFont(62)
uB5Hist.GetYaxis().SetTitleOffset(0.27)
uB5Hist.Draw("C")
lB5Hist.Draw("SAMEC")
tex.DrawLatexNDC(0.05,0.15,"b_{2}^{(1)}")

canComp.cd()
padB7 = ROOT.TPad("padB7", "padB7", BxPad, 1 - 4*ypixRat - ByPad, Bx1E, 1 - 3*ypixRat - ByPad)
padB7.SetRightMargin(0.0)
padB7.SetLeftMargin(0.0)
padB7.SetTopMargin(0.0)
padB7.SetBottomMargin(0.0)
padB7.SetFrameFillColor(0)
padB7.SetBorderSize(0)
padB7.Draw()
padB7.cd()
uB6Hist.GetXaxis().SetLabelSize(labelSize*3)
uB6Hist.GetXaxis().SetTitleSize(labelSize*3)
uB6Hist.GetXaxis().SetLabelFont(62)
uB6Hist.GetXaxis().SetTitleFont(62)
uB6Hist.GetXaxis().SetTitleOffset(0.91)
uB6Hist.GetXaxis().SetLabelColor(uColor)
uB6Hist.GetYaxis().SetLabelSize(labelSize*3)
uB6Hist.GetYaxis().SetTitleSize(labelSize*3)
uB6Hist.GetYaxis().SetLabelFont(62)
uB6Hist.GetYaxis().SetTitleFont(62)
uB6Hist.GetYaxis().SetTitleOffset(0.27)
uB6Hist.Draw("C")
lB6Hist.Draw("SAMEC")
tex.DrawLatexNDC(0.05,0.15,"b_{3}^{(0)}")

canComp.cd()
padB8 = ROOT.TPad("padB8", "padB8", Bx2S, 1 - 4*ypixRat - ByPad, xpixRat, 1 -3*ypixRat - ByPad)
padB8.SetRightMargin(0.0)
padB8.SetLeftMargin(0.0)
padB8.SetTopMargin(0.0)
padB8.SetBottomMargin(0.0)
padB8.SetFrameFillColor(0)
padB8.SetBorderSize(0)
padB8.Draw()
padB8.cd()
uB7Hist.GetXaxis().SetLabelSize(labelSize*3)
uB7Hist.GetXaxis().SetTitleSize(labelSize*3)
uB7Hist.GetXaxis().SetLabelFont(62)
uB7Hist.GetXaxis().SetTitleFont(62)
uB7Hist.GetXaxis().SetTitleOffset(0.91)
uB7Hist.GetXaxis().SetLabelColor(uColor)
uB7Hist.GetYaxis().SetLabelSize(0)
uB7Hist.GetYaxis().SetTitleSize(0)
uB7Hist.GetYaxis().SetLabelFont(62)
uB7Hist.GetYaxis().SetTitleFont(62)
uB7Hist.GetYaxis().SetTitleOffset(0.27)
uB7Hist.Draw("C")
lB7Hist.Draw("SAMEC")
tex.DrawLatexNDC(0.05,0.15,"b_{3}^{(1)}")
"""

"""
bsXPos = 0.13
bsXLen = 0.17
bsYLen = 0.25
yShift = 0.01
xShift = 0.044
canComp.cd()
blank = ROOT.TPad("padB2", "padB2", 0, 0 + yShift, xShift, 1 - ypixRat + yShift)
blank.SetFrameFillColor(0)
blank.SetBorderSize(0)
blank.Draw()
"""


print("here2")
canComp.cd()
canComp.Update()
canComp.Print("compareTimeBasis.png")
print("here3")






#####################
##  Compare Coeff  ##
#####################

uCoefImg = np.zeros((Npix, Npix), dtype = np.double)
lCoefImg = np.zeros((Npix, Npix), dtype = np.double)
#for ilg in range(1,Nlg):
for ilg in range(Nlg):
  lCoefImg[:,:] += lfftCoef[ilg,0]*legImg[ilg,:,:]
  if ilg is not 0:
    uCoefImg[:,:] += ufftCoef[ilg,0,:,:]*legImg[ilg,:,:]

uNorm = 1.0/max(np.absolute([np.amax(uCoefImg), np.amin(uCoefImg)]))
lNorm = 1.0/max(np.absolute([np.amax(lCoefImg), np.amin(lCoefImg)]))

uCoefImg *= uNorm
lCoefImg *= lNorm

uCoefHist = ROOT.TH2F("uCoefHist", "uCoefHist", Npix/2+1, qmin, qmax, Npix, qmin, qmax)
lCoefHist = ROOT.TH2F("lCoefHist", "lCoefHist", Npix/2+1, qmin, qmax, Npix, qmin, qmax)
uCoefHist.SetContour(80)
lCoefHist.SetContour(80)

for ir in range(Npix):
  for ic in range(Npix/2 + 1):
    lCoefHist.SetBinContent(ic+1, ir+1, lCoefImg[ir,ic])
    uCoefHist.SetBinContent(ic+1, ir+1, uCoefImg[ir,Npix/2 + ic])

for ir in range(Npix):
  for ic in range(Npix/2 + 1):
    if np.sqrt((ir-Npix/2)**2 + (ic-Npix/2)**2) > Npix/2:
      uCoefHist.SetBinContent(Npix/2 - ic + 2, ir+1, -1e20)
      lCoefHist.SetBinContent(ic+1, ir+1, -1e20)

uCoefHist.SetMaximum(1)
lCoefHist.SetMaximum(1)
uCoefHist.SetMinimum(-1)
lCoefHist.SetMinimum(-1)

uCoefHist.GetXaxis().SetLabelSize(0)
uCoefHist.GetYaxis().SetLabelSize(0)
uCoefHist.GetXaxis().SetTickSize(0)
uCoefHist.GetYaxis().SetTickSize(0)
uCoefHist.GetYaxis().SetTitleColor(0)
lCoefHist.GetXaxis().SetLabelSize(0)
lCoefHist.GetYaxis().SetLabelSize(0)
lCoefHist.GetXaxis().SetTickSize(0)
lCoefHist.GetYaxis().SetTickSize(0)
lCoefHist.GetYaxis().SetTitleColor(0)




xpix = 900.
ypix = 1000.
ycut = 270.
canCoeff = ROOT.TCanvas("canCoeff", "canCoeff", int(xpix), int(ypix))
canCoeff.cd()

pOff = 0.03
pSep = 0.01
pSize = 0.4
canCoeff.cd()
print(ycut/ypix, (ypix-ycut)/xpix/2)
ccpad1 = ROOT.TPad("ccpad1", "ccpad1", pOff, ycut/ypix, pOff + (ypix-ycut)/xpix/2, 1)
ccpad1.SetRightMargin(0)
ccpad1.SetLeftMargin(0)
ccpad1.SetTopMargin(0)
ccpad1.SetBottomMargin(0)
ccpad1.SetFrameFillColor(0)
ccpad1.SetFrameLineColor(0)
ccpad1.SetBorderSize(0)
ccpad1.Draw()
ccpad1.cd()
lCoefHist.Draw("COL")
lcirc.Draw("SAME")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)
lAxis.Draw("SAME")
txt.SetTextSize(labelSize*1.5)
txt.SetTextFont(62)
txt.SetTextColor(1)

for il in range(len(laxisLabels)):
  txt.DrawText(il*laxisPosItr, 0.447, laxisLabels[il])
txt.SetTextSize(labelSize*1.5)
txt.DrawText(0.041, 0.41, "E [eV]")



thist = ROOT.TH2F("thist", "thist", 1, 0, 1, 1, 0, 1)
thist.SetBinContent(1,1,-5)
thist.SetMaximum(1)
thist.SetMinimum(-1)
thist.SetContour(80)
thist.GetXaxis().SetLabelSize(0)
thist.GetXaxis().SetTickSize(0)
thist.GetYaxis().SetTitleColor(0)
thist.GetYaxis().SetLabelSize(0)
thist.GetYaxis().SetTickSize(0)

canCoeff.cd()
ccpadFL = ROOT.TPad("ccpadFL", "ccpadFL", (ypix-ycut)/xpix/2 + pSep + pOff, ycut/ypix, 1, 1)
ccpadFL.SetRightMargin((1 - (ypix-ycut)/xpix + pSep)/(1 - ((ypix-ycut)/xpix/2 + pSep))*0.55)
ccpadFL.SetLeftMargin(0)
ccpadFL.SetTopMargin(0)
ccpadFL.SetBottomMargin(0)
ccpadFL.SetFrameFillColor(0)
ccpadFL.SetFrameLineColor(0)
ccpadFL.SetBorderSize(0)
ccpadFL.Draw()
ccpadFL.cd()
thist.Draw("COLZ")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC((1 - (1 - (ypix-ycut)/xpix + pSep)/(1 - ((ypix-ycut)/xpix/2 + pSep))*0.55), 0, (1 - 0.55*(1 - (ypix-ycut)/xpix + pSep)/(1 - ((ypix-ycut)/xpix/2 + pSep))), 1)
tb.DrawLineNDC(0, 1, 1 - 0.55*(1 - (ypix-ycut)/xpix + pSep)/(1 - ((ypix-ycut)/xpix/2 + pSep)), 1)
bb.DrawLineNDC(0, 0, 1 - 0.55*(1 - (ypix-ycut)/xpix + pSep)/(1 - ((ypix-ycut)/xpix/2 + pSep)), 0)

canCoeff.cd()
ccpad2 = ROOT.TPad("ccpad2", "ccpad2", (ypix-ycut)/xpix/2 + pSep + pOff, ycut/ypix, (ypix-ycut)/xpix + pSep + pOff, 1)
ccpad2.SetRightMargin(0)
ccpad2.SetLeftMargin(0)
ccpad2.SetTopMargin(0)
ccpad2.SetBottomMargin(0)
ccpad2.SetFrameFillColor(0)
ccpad2.SetFrameLineColor(0)
ccpad2.SetBorderSize(0)
ccpad2.Draw()
ccpad2.cd()
uCoefHist.Draw("COL")
ucirc.Draw("SAME")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)
uAxis.Draw("SAME")





######################
##  Power Spectrum  ##
######################

NbasisPlot = 8
uPowSpect = np.zeros((NfftBasis,1), dtype=np.double)
lPowSpect = np.zeros((NfftBasis,1), dtype=np.double)
#for ilg in range(1, Nlg):
for ilg in range(Nlg):
  #if ilg >2:
  #  continue
  #uWeight = np.sum(ufftCoefRad[ilg, :, :]**2)
  #lWeight = np.sum(lfftCoefRad[ilg, :, :]**2)
  for ibs in range(NbasisPlot*2):
    for i in range(NuRad):
      #if i < 20:
      #  continue
      if ilg is not 0:
        uPowSpect[ibs] += ufftCoefRad[ilg, ibs, i]**2
      if ilg == 2:
        print("testing vals",ibs,i,ufftCoefRad[ilg, ibs, i]**2)
    for i in range(10,15):
      lPowSpect[ibs] += lfftCoefRad[ilg, ibs, i]**2
print(uPowSpect)
print(sum(uPowSpect))
uPowSpect /= sum(uPowSpect)
lPowSpect /= sum(lPowSpect)
for i in range(NbasisPlot):
  print(lPowSpect[i],uPowSpect[i])

uPowHistF = ROOT.TH1F("uPowHistF", "uPowHistF", NbasisPlot, 0, NbasisPlot)
uPowHistD = ROOT.TH1F("uPowHistD", "uPowHistF", NbasisPlot, 0, NbasisPlot)
lPowHistF = ROOT.TH1F("lPowHistF", "pPowHistF", NbasisPlot, 0, NbasisPlot)
lPowHistD = ROOT.TH1F("lPowHistD", "pPowHistD", NbasisPlot, 0, NbasisPlot)
for i in range(NbasisPlot):
  uPowHistF.SetBinContent(i+1, uPowSpect[2*i])
  lPowHistF.SetBinContent(i+1, lPowSpect[2*i])
  uPowHistD.SetBinContent(i+1, uPowSpect[2*i+1])
  lPowHistD.SetBinContent(i+1, lPowSpect[2*i+1])
  """
  uPowHistF.SetBinContent(i+1, uPowSpect[i])
  lPowHistF.SetBinContent(i+1, lPowSpect[i])
  uPowHistD.SetBinContent(i+1, uPowSpect[i+NfftPows])
  lPowHistD.SetBinContent(i+1, lPowSpect[i+NfftPows])
  """
  uPowHistF.SetBinError(i+1, 1e-10)
  lPowHistF.SetBinError(i+1, 1e-10)
  uPowHistD.SetBinError(i+1, 1e-10)
  lPowHistD.SetBinError(i+1, 1e-10)
  print("pows",i,uPowSpect[2*i], uPowSpect[2*i+1], lPowSpect[2*i], lPowSpect[2*i+1])

uPowHistF.GetXaxis().SetTitle("Basis Index")
uPowHistF.GetYaxis().SetTitle("Power")
uPowHistF.GetYaxis().SetTitleSize(0.0)
uPowHistF.GetXaxis().SetTitleSize(0.1)
uPowHistF.GetXaxis().SetLabelSize(0.08)
uPowHistF.GetYaxis().SetLabelSize(0.0)
uPowHistF.GetYaxis().SetTitleOffset(0.33)
uPowHistF.GetXaxis().SetTitleOffset(0.9)
uPowHistF.GetXaxis().SetLabelFont(62)
uPowHistF.GetYaxis().SetLabelFont(62)
uPowHistF.GetXaxis().SetTitleFont(62)
uPowHistF.GetYaxis().CenterTitle(1)
uPowHistF.GetXaxis().CenterTitle(1)
uPowHistF.GetXaxis().CenterLabels(1)
uPowHistF.SetTitle("UED Time Basis Power Spectrum")

#uPowHistF.SetMarkerColor(uColor)
uPowHistF.SetMarkerSize(0)
uPowHistF.SetLineColor(uColor)
uPowHistD.SetMarkerSize(0)
uPowHistD.SetLineColor(uColor)
uPowHistD.SetLineStyle(7)

lPowHistF.GetXaxis().SetTitle("Basis Index")
lPowHistF.GetYaxis().SetTitleSize(0)
lPowHistF.GetXaxis().SetTitleSize(0.1)
lPowHistF.GetXaxis().SetLabelSize(0.08)
lPowHistF.GetYaxis().SetLabelSize(0.08)
lPowHistF.GetYaxis().SetTitleOffset(0.33)
lPowHistF.GetXaxis().SetTitleOffset(0.9)
lPowHistF.GetXaxis().SetLabelFont(62)
lPowHistF.GetYaxis().SetLabelFont(62)
lPowHistF.GetXaxis().SetTitleFont(62)
lPowHistF.GetYaxis().CenterTitle(1)
lPowHistF.GetXaxis().CenterTitle(1)
lPowHistF.GetXaxis().CenterLabels(1)
#lPowHistF.SetMarkerColor(uColor)
lPowHistF.SetMarkerSize(0)
lPowHistF.SetLineColor(lColor)
lPowHistD.SetMarkerSize(0)
lPowHistD.SetLineColor(lColor)
lPowHistD.SetLineStyle(7)

# Setting same minimum
minimum = 0.9*min([lPowHistF.GetMinimum(), lPowHistD.GetMinimum(),
                uPowHistF.GetMinimum(), uPowHistD.GetMinimum()])

uPowHistF.SetMaximum(1);
uPowHistD.SetMaximum(1);
lPowHistF.SetMaximum(1);
lPowHistD.SetMaximum(1);
uPowHistF.SetMinimum(minimum);
uPowHistD.SetMinimum(minimum);
lPowHistF.SetMinimum(minimum);
lPowHistD.SetMinimum(minimum);


lPowLgnd = ROOT.TLegend(0.72, 0.6, 1.18, 0.9)
lPowLgnd.SetFillStyle(0)
lPowLgnd.SetBorderSize(0)
lPowLgnd.SetTextSize(0.09)
lPowLgnd.AddEntry(lPowHistF, "b_{index}", "l")
lPowLgnd.AddEntry(lPowHistD, "b_{index}^{(d)}", "l")

uPowLgnd = ROOT.TLegend(0.5, 0.6, 0.88, 0.9)
uPowLgnd.SetFillStyle(0)
uPowLgnd.SetBorderSize(0)
uPowLgnd.SetTextSize(0.09)
uPowLgnd.AddEntry(uPowHistF, "b_{index}", "l")
uPowLgnd.AddEntry(uPowHistD, "b_{index}^{(d)}", "l")


canCoeff.cd()
ccpad3 = ROOT.TPad("ccpad3", "ccpad3", 0, 0, pOff + (ypix-ycut)/xpix/2, (ycut-20)/ypix)
ccpad3.SetRightMargin(0.005)
ccpad3.SetLeftMargin(0.09)
ccpad3.SetTopMargin(0.1)
ccpad3.SetBottomMargin(0.2)
ccpad3.SetFrameFillColor(0)
ccpad3.SetFrameLineColor(0)
ccpad3.SetBorderSize(0)
ccpad3.SetLogy(1)
ccpad3.Draw()
ccpad3.cd()

lPowHistF.SetMaximum(1)
lPowHistF.Draw("LPE")
lPowHistD.Draw("SAMELPE")
lPowLgnd.Draw("SAME")

txt.SetTextSize(labelSize*3)
txt.DrawText(4, 1.7, "LCLS Power Spectrum")

canCoeff.cd()
ccpad4 = ROOT.TPad("ccpad4", "ccpad4", (ypix-ycut)/xpix/2 + pOff + pSep, 0, 1, (ycut-20)/ypix)
ccpad4.SetRightMargin(0.09 + (1 - (ypix-ycut)/xpix + pSep)/(1 - ((ypix-ycut)/xpix/2 + pSep))*0.55)
ccpad4.SetLeftMargin(0.005)
ccpad4.SetTopMargin(0.1)
ccpad4.SetBottomMargin(0.2)
ccpad4.SetFrameFillColor(0)
ccpad4.SetFrameLineColor(0)
ccpad4.SetBorderSize(0)
ccpad4.SetLogy(1)
ccpad4.Draw()
ccpad4.cd()

uPowHistF.SetMaximum(1)
uPowHistF.Draw("LPE")
uPowHistD.Draw("SAMEE")
uPowLgnd.Draw("SAME")

txt.DrawText(4, 1.7, "UED Power Spectrum")

tex = ROOT.TLatex()
tex.SetTextSize(0.12)
tex.SetTextFont(42)
axLabels = ["b_{0}^{(0)}", "b_{0}^{(1)}", "b_{1}^{(0)}", "b_{1}^{(1)}", "b_{2}^{(0)}", "b_{2}^{(1)}"]
#axLabels = ["b_{0}^{(0)}", "b_{0}^{(1)}", "b_{1}^{(0)}", "b_{1}^{(1)}", "b_{2}^{(0)}", "b_{2}^{(1)}", "b_{3}^{(0)}", "b_{3}^{(1)}", "b_{4}^{(0)}", "b_{4}^{(1)}", "b_{5}^{(0)}", "b_{5}^{(1)}"]

#print(uPowHist.GetMinimum())
#for i in range(len(axLabels)):
  #tex.DrawLatexNDC(0.115+0.0683*i, 0.155, axLabels[i])
  #tex.DrawLatexNDC(0.145+0.0683*2*i, 0.155, axLabels[i])

canCoeff.Print("compareCoeff.png")




##################################
#####  Comparing Time Basis  #####
##################################

uedOrthoBasis = np.fromfile("../rotorBasis/uedfftRotorbasisFine.dat", dtype=np.double)
uedSemiOrthoBasis = np.fromfile("../rotorBasis/uedfftRotorbasisFine_semiOrthog.dat", dtype=np.double)
lclsOrthoBasis = np.fromfile("../rotorBasis/lclsfftRotorbasis.dat", dtype=np.double)
lclsSemiOrthoBasis = np.fromfile("../rotorBasis/lclsfftRotorbasis_semiOrthog.dat", dtype=np.double)

print("Northo basis",2*NfftBasis)
uedOrthoBasis = np.reshape(uedOrthoBasis, (NfftBasis, -1))
uedSemiOrthoBasis = np.reshape(uedSemiOrthoBasis, (NfftBasis, -1))
lclsOrthoBasis = np.reshape(lclsOrthoBasis, (NfftBasis, -1))
lclsSemiOrthoBasis = np.reshape(lclsSemiOrthoBasis, (NfftBasis, -1))

"""
uedOrthoBasis *= 3
uedSemiOrthoBasis *= 3
lclsOrthoBasis *= 3
lclsSemiOrthoBasis *= 3
"""

ylabels = []
yticks = []
xticks = [317, 318, 319, 320]
uedX = np.linspace(315.914, 320.214, uedOrthoBasis.shape[1])
lclsX = np.linspace(315.914, 320.214, lclsOrthoBasis.shape[1])
fig, ax = plt.subplots(1,4, figsize=(16,7))
for i in range(min([5,NfftBasis/2])):
  ax[0].plot(uedX, uedSemiOrthoBasis[2*i,:] + 2*i, "-b")
  ax[0].plot(uedX, uedSemiOrthoBasis[2*i+1,:] + (2*i+1), "b--")
  ax[1].plot(uedX, uedOrthoBasis[2*i,:] + 2*i, "-b")
  ax[1].plot(uedX, uedOrthoBasis[2*i+1,:] + (2*i+1), "b--")
  ax[2].plot(lclsX, lclsSemiOrthoBasis[2*i,:] + 2*i, "-r")
  ax[2].plot(lclsX, lclsSemiOrthoBasis[2*i+1,:] + (2*i+1), "r--")
  ax[3].plot(lclsX, lclsOrthoBasis[2*i,:] + 2*i, "-r")
  ax[3].plot(lclsX, lclsOrthoBasis[2*i+1,:] + (2*i+1), "r--")
  yticks.append(2*i)
  yticks.append(2*i+1)
  yl1 = "$\mathregular{b_"+str(i)+"}$"
  yl2 = "$\mathregular{b^{d}_"+str(i)+"}$"
  ylabels.append(yl1)
  ylabels.append(yl2)
ax[0].set_xlim([316.5, 320])
ax[1].set_xlim([316.5, 320])
ax[2].set_xlim([316.5, 320])
ax[3].set_xlim([316.5, 320])
ax[0].set_ylim([-0.5,2*(min([5,NfftBasis/2])-1)+1.5])
ax[1].set_ylim([-0.5,2*(min([5,NfftBasis/2])-1)+1.5])
ax[2].set_ylim([-0.5,2*(min([5,NfftBasis/2])-1)+1.5])
ax[3].set_ylim([-0.5,2*(min([5,NfftBasis/2])-1)+1.5])

ax[0].tick_params(labelsize=14)
ax[1].tick_params(labelsize=14)
ax[2].tick_params(labelsize=14)
ax[3].tick_params(labelsize=14)
ax[0].set_xticks(xticks)
ax[1].set_xticks(xticks)
ax[2].set_xticks(xticks)
ax[3].set_xticks(xticks)
ax[0].set_yticks(yticks)
ax[0].set_yticklabels(ylabels)
ax[1].set_yticklabels([])
ax[2].set_yticklabels([])
ax[3].set_yticklabels([])

ax[0].set_xlabel("Time [ps]", fontsize=17)
ax[1].set_xlabel("Time [ps]", fontsize=17)
ax[2].set_xlabel("Time [ps]", fontsize=17)
ax[3].set_xlabel("Time [ps]", fontsize=17)

ax[0].text(0.91, 0.985, "a", color='k', 
    transform=ax[0].transAxes, verticalalignment='top', 
    fontsize=27)
ax[1].text(0.91, 0.985, "b", color='k', 
    transform=ax[1].transAxes, verticalalignment='top', 
    fontsize=27)
ax[2].text(0.91, 0.985, "c", color='k', 
    transform=ax[2].transAxes, verticalalignment='top', 
    fontsize=27)
ax[3].text(0.91, 0.985, "d", color='k', 
    transform=ax[3].transAxes, verticalalignment='top', 
    fontsize=27)
fig.tight_layout()
plt.savefig("timeBasisSet.png")
plt.close()




"""
pnc bank
channel services
500 first ave
pittsburg pa, 15219
mailstop p7-pfsc-01-s
"""

