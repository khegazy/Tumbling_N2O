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

maxQ = 14.1207
circRad = 0.14
labelSize = 0.04

####  Canvas ####
xpix = 820.
ypix = 700.
can = ROOT.TCanvas("can", "can", int(xpix), int(ypix) + 28)
can.cd()


dataFile = ROOT.TFile.Open("../UED/output/data/combinedResults.root", "read")
simFile = ROOT.TFile.Open("../UED/output/MC/combinedResults.root", "read")




dhist = dataFile.Get("diffSMS26").Clone("dclone")
shist = simFile.Get("diffSMS26").Clone("dclone")
Npix = dhist.GetNbinsX()
maxQ = dhist.GetXaxis().GetXmax()



dDhist = ROOT.TH2F("dDhist", "dDhist", Npix/2+1, 0, 1, Npix, 0, 1)
sDhist = ROOT.TH2F("sDhist", "sDhist", Npix/2+1, 0, 1, Npix, 0, 1)
dDhist.SetContour(80)
sDhist.SetContour(80)

for ir in range(Npix):
  for ic in range(Npix/2 + 1):
    dDhist.SetBinContent(ic+1, ir+1, dhist.GetBinContent(ic+1, ir+1))
    sDhist.SetBinContent(ic+1, ir+1, shist.GetBinContent(Npix/2 + ic, ir+1))
    if (ir - Npix/2)**2 + (Npix/2 - ic)**2 < 400:
      dDhist.SetBinContent(ic+1, ir+1, 0)
    if (ir - Npix/2)**2 + (ic)**2 < 400:
      sDhist.SetBinContent(ic+1, ir+1, 0)


dDhist.Scale(5/np.max(np.abs([dDhist.GetMaximum(), dDhist.GetMinimum()])))
sDhist.Scale(1./np.max(np.abs([sDhist.GetMaximum(), sDhist.GetMinimum()])))

for ir in range(Npix):
  for ic in range(Npix/2 + 1):
    if np.sqrt((ir-Npix/2)**2 + (ic-Npix/2)**2) > Npix/2:
      sDhist.SetBinContent(Npix/2 - ic + 2, ir+1, -1e20)
      dDhist.SetBinContent(ic+1, ir+1, -1e20)

dDhist.SetMaximum(1)
sDhist.SetMaximum(1)
dDhist.SetMinimum(-1)
sDhist.SetMinimum(-1)

dDhist.GetXaxis().SetLabelSize(0)
dDhist.GetYaxis().SetLabelSize(0)
dDhist.GetXaxis().SetTickSize(0)
dDhist.GetYaxis().SetTickSize(0)
dDhist.GetYaxis().SetTitleColor(0)
sDhist.GetXaxis().SetLabelSize(0)
sDhist.GetYaxis().SetLabelSize(0)
sDhist.GetXaxis().SetTickSize(0)
sDhist.GetYaxis().SetTickSize(0)
sDhist.GetYaxis().SetTitleColor(0)

lcirc = ROOT.TEllipse(1.0, 0.5, circRad, circRad/2, 90, 270)
lcirc.SetFillColor(1)
rcirc = ROOT.TEllipse(0.0, 0.5, circRad, circRad/2)
rcirc.SetFillColor(1)

rAxis = ROOT.TGaxis(circRad, 0.5, 1.0, 0.5, maxQ*circRad, maxQ, 507, "-+")
rAxis.SetLabelSize(labelSize*1.5)
rAxis.SetLabelOffset(-0.008)
rAxis.SetTitleSize(labelSize*1.5)
rAxis.SetTitle("Q [#AA^{-1}]")
lAxis = ROOT.TGaxis(0, 0.5, 1 - circRad, 0.5, maxQ, maxQ*circRad, 507, "-+")
lAxis.SetLabelSize(0)
laxisLabels = ["14", "12", "10", "8", "6", "4", "2"]
laxisPosItr = (1.0-circRad)/(len(laxisLabels) - 1)

lb = ROOT.TLine()
rb = ROOT.TLine()
tb = ROOT.TLine()
bb = ROOT.TLine()
lb.SetLineColor(0)
rb.SetLineColor(0)
tb.SetLineColor(0)
bb.SetLineColor(0)

txt = ROOT.TText()
txt.SetTextSize(labelSize)
txt.SetTextAlign(22)
txt.SetTextColor(0)

tex = ROOT.TLatex()
tex.SetTextFont(62)
tex.SetTextSize(labelSize)


pOff = 0.03
pSep = 0.01
pSize = 0.4
can.cd()
pad1 = ROOT.TPad("pad1", "pad1", pOff, 0.02, pOff + ypix/xpix/2, 1 - 0.02)
pad1.SetRightMargin(0)
pad1.SetLeftMargin(0)
pad1.SetTopMargin(0)
pad1.SetBottomMargin(0)
pad1.SetFrameFillColor(0)
pad1.SetFrameLineColor(0)
pad1.SetBorderSize(0)
pad1.Draw()
pad1.cd()
dDhist.Draw("COL")
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
  txt.DrawText(il*laxisPosItr, 0.445, laxisLabels[il])
tex.SetTextSize(labelSize*1.5)
tex.DrawLatexNDC(0.0, 0.395, "Q [#AA^{-1}]")
txt.SetTextSize(labelSize*5)
txt.DrawTextNDC(0.1, 0.9, "a")


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

can.cd()
pad2a = ROOT.TPad("pad2a", "pad2a", ypix/xpix/2 + pSep + pOff, 0.02, 1, 1 - 0.02)
pad2a.SetRightMargin((1 - ypix/xpix + pSep)/(1 - (ypix/xpix/2 + pSep))*0.55)
pad2a.SetLeftMargin(0)
pad2a.SetTopMargin(0)
pad2a.SetBottomMargin(0)
pad2a.SetFrameFillColor(0)
pad2a.SetFrameLineColor(0)
pad2a.SetBorderSize(0)
pad2a.Draw()
pad2a.cd()
thist.Draw("COLZ")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC((1 - (1 - ypix/xpix + pSep)/(1 - (ypix/xpix/2 + pSep))*0.55), 0, (1 - 0.55*(1 - ypix/xpix + pSep)/(1 - (ypix/xpix/2 + pSep))), 1)
tb.DrawLineNDC(0, 1, 1 - 0.55*(1 - ypix/xpix + pSep)/(1 - (ypix/xpix/2 + pSep)), 1)
bb.DrawLineNDC(0, 0, 1 - 0.55*(1 - ypix/xpix + pSep)/(1 - (ypix/xpix/2 + pSep)), 0)

can.cd()
pad2 = ROOT.TPad("pad2", "pad2", ypix/xpix/2 + pSep + pOff, 0.02, ypix/xpix + pSep + pOff, 1 - 0.02)
pad2.SetRightMargin(0)
pad2.SetLeftMargin(0)
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0)
pad2.SetFrameFillColor(0)
pad2.SetFrameLineColor(0)
pad2.SetBorderSize(0)
pad2.Draw()
pad2.cd()
sDhist.Draw("COL")
rcirc.Draw("SAME")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)
rAxis.Draw("SAME")

txt.SetTextSize(labelSize*5)
txt.DrawTextNDC(0.9, 0.9, "b")

can.Print("diffPmaxCompare.png")



####  Canvas ####
canS = ROOT.TCanvas("canS", "canS", 850, 850)
canS.cd()

canS.SetRightMargin(0)
canS.SetLeftMargin(0)
canS.SetTopMargin(0)
canS.SetBottomMargin(0)
canS.SetFrameFillColor(0)
canS.SetFrameLineColor(0)
canS.SetBorderSize(0)


dPatt = ROOT.TH2F("dPatt", "dPatt", Npix, 0, 1, Npix, 0, 1)
dPatt.SetContour(80)

for ir in range(Npix):
  for ic in range(Npix):
    dPatt.SetBinContent(ic+1, ir+1, dhist.GetBinContent(ic+1, ir+1))
    if (ir - Npix/2)**2 + (Npix/2 - ic)**2 < 400:
      dPatt.SetBinContent(ic+1, ir+1, 0)

dPatt.Scale(1./np.max(np.abs([dPatt.GetMaximum(), dPatt.GetMinimum()])))

for ir in range(Npix):
  for ic in range(Npix):
    if np.sqrt((ir-Npix/2)**2 + (ic-Npix/2)**2) > Npix/2:
      dPatt.SetBinContent(ic+1, ir+1, -1e20)


dPatt.SetMaximum(1)
dPatt.SetMinimum(-1)

dPatt.GetXaxis().SetLabelSize(0)
dPatt.GetYaxis().SetLabelSize(0)
dPatt.GetXaxis().SetTickSize(0)
dPatt.GetYaxis().SetTickSize(0)
dPatt.GetYaxis().SetTitleColor(0)

lcirc = ROOT.TEllipse(0.5, 0.5, circRad/2., circRad/2.)
lcirc.SetFillColor(1)

dPatt.Draw("COL")
lcirc.Draw("SAME")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)

canS.Print("setupPatt.png")





######################
##  AutoCorr Plots  ##
######################

canA = ROOT.TCanvas("canA", "canA", int(xpix), int(ypix) + 28)
canA.cd()

dhist = dataFile.Get("autCRel26").Clone("dclone")
shist = simFile.Get("autCRel26").Clone("dclone")
Npix = dhist.GetNbinsX()
maxAA = dhist.GetXaxis().GetXmax()
print(maxAA)

dAChist = ROOT.TH2F("dAChist", "dAChist", Npix/2+1, 0, maxAA, Npix, -maxAA, maxAA)
sAChist = ROOT.TH2F("sAChist", "sAChist", Npix/2+1, 0, maxAA, Npix, -maxAA, maxAA)
dAChist.SetContour(80)
sAChist.SetContour(80)

for ir in range(Npix):
  for ic in range(Npix/2 + 1):
    dAChist.SetBinContent(ic+1, ir+1, dhist.GetBinContent(ic+1, ir+1))
    sAChist.SetBinContent(ic+1, ir+1, shist.GetBinContent(Npix/2 + ic, ir+1))
    if (ir - Npix/2)**2 + (Npix/2 - ic)**2 < 400:
      dAChist.SetBinContent(ic+1, ir+1, 0)
    if (ir - Npix/2)**2 + (ic)**2 < 400:
      sAChist.SetBinContent(ic+1, ir+1, 0)


dAChist.Scale(1./np.max(np.abs([dAChist.GetMaximum(), dAChist.GetMinimum()])))
sAChist.Scale(1./np.max(np.abs([sAChist.GetMaximum(), sAChist.GetMinimum()])))

for ir in range(Npix):
  for ic in range(Npix/2 + 1):
    if np.sqrt((ir-Npix/2)**2 + (ic-Npix/2)**2) > Npix/2:
      sAChist.SetBinContent(Npix/2 - ic + 2, ir+1, -1e20)
      dAChist.SetBinContent(ic+1, ir+1, -1e20)

dAChist.SetMaximum(1)
sAChist.SetMaximum(1)
dAChist.SetMinimum(-1)
sAChist.SetMinimum(-1)

dAChist.GetXaxis().SetLabelSize(0)
dAChist.GetYaxis().SetLabelSize(0)
dAChist.GetXaxis().SetTickSize(0)
dAChist.GetYaxis().SetTickSize(0)
dAChist.GetYaxis().SetTitleColor(0)
sAChist.GetXaxis().SetLabelSize(0)
sAChist.GetYaxis().SetLabelSize(0)
sAChist.GetXaxis().SetTickSize(0)
sAChist.GetYaxis().SetTickSize(0)
sAChist.GetYaxis().SetTitleColor(0)

lcirc = ROOT.TEllipse(maxAA, 0, 0.5, 0.5, 90, 270)
lcirc.SetFillColor(1)
rcirc = ROOT.TEllipse(0.0, 0, 0.5, 0.5)
rcirc.SetFillColor(1)

rAxis = ROOT.TGaxis(0.5, 0, maxAA, 0, 0.5, maxAA, 507, "-+")
rAxis.SetLabelSize(labelSize*1.5)
rAxis.SetLabelOffset(-0.008)
rAxis.SetTitleSize(labelSize*1.5)
rAxis.SetTitle("[#AA]")
lShift = 0.18
lAxis = ROOT.TGaxis(0, 0, maxAA - 0.5, 0, maxAA, 0.5, 507, "-+")
lAxis.SetLabelSize(0)
laxisLabels = ["2.5", "2", "1.5", "1", "0.5"]
laxisPosItr = ((maxAA-lShift)-0.5)/(len(laxisLabels) - 1)

lb = ROOT.TLine()
rb = ROOT.TLine()
tb = ROOT.TLine()
bb = ROOT.TLine()
lb.SetLineColor(0)
rb.SetLineColor(0)
tb.SetLineColor(0)
bb.SetLineColor(0)

txt = ROOT.TText()
txt.SetTextSize(labelSize)
txt.SetTextAlign(22)
txt.SetTextColor(0)



pOff = 0.03
pSep = 0.01
pSize = 0.4
canA.cd()
pad1 = ROOT.TPad("pad1", "pad1", pOff, 0.02, pOff + ypix/xpix/2, 1 - 0.02)
pad1.SetRightMargin(0)
pad1.SetLeftMargin(0)
pad1.SetTopMargin(0)
pad1.SetBottomMargin(0)
pad1.SetFrameFillColor(0)
pad1.SetFrameLineColor(0)
pad1.SetBorderSize(0)
pad1.Draw()
pad1.cd()
dAChist.Draw("COL")
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
  txt.DrawText(lShift + il*laxisPosItr, -2*(0.5-0.445)*maxAA, laxisLabels[il])
tex.SetTextSize(labelSize*1.5)
tex.DrawLatexNDC(0.0,  0.395, "[#AA]")
txt.SetTextSize(labelSize*5)
txt.DrawTextNDC(0.1, 0.9, "c")



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

canA.cd()
pad2a = ROOT.TPad("pad2a", "pad2a", ypix/xpix/2 + pSep + pOff, 0.02, 1, 1 - 0.02)
pad2a.SetRightMargin((1 - ypix/xpix + pSep)/(1 - (ypix/xpix/2 + pSep))*0.55)
pad2a.SetLeftMargin(0)
pad2a.SetTopMargin(0)
pad2a.SetBottomMargin(0)
pad2a.SetFrameFillColor(0)
pad2a.SetFrameLineColor(0)
pad2a.SetBorderSize(0)
pad2a.Draw()
pad2a.cd()
thist.Draw("COLZ")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC((1 - (1 - ypix/xpix + pSep)/(1 - (ypix/xpix/2 + pSep))*0.55), 0, (1 - 0.55*(1 - ypix/xpix + pSep)/(1 - (ypix/xpix/2 + pSep))), 1)
tb.DrawLineNDC(0, 1, 1 - 0.55*(1 - ypix/xpix + pSep)/(1 - (ypix/xpix/2 + pSep)), 1)
bb.DrawLineNDC(0, 0, 1 - 0.55*(1 - ypix/xpix + pSep)/(1 - (ypix/xpix/2 + pSep)), 0)

canA.cd()
pad2 = ROOT.TPad("pad2", "pad2", ypix/xpix/2 + pSep + pOff, 0.02, ypix/xpix + pSep + pOff, 1 - 0.02)
pad2.SetRightMargin(0)
pad2.SetLeftMargin(0)
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0)
pad2.SetFrameFillColor(0)
pad2.SetFrameLineColor(0)
pad2.SetBorderSize(0)
pad2.Draw()
pad2.cd()
sAChist.Draw("COL")
rcirc.Draw("SAME")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)
rAxis.Draw("SAME")
txt.SetTextSize(labelSize*5)
txt.DrawTextNDC(0.9, 0.9, "d")


canA.Print("autCmaxCompare.png")




































"""

xpix = 1500
ypix = 600
canComp = ROOT.TCanvas("canComp", "canComp", xpix, ypix);

canComp.SetRightMargin(0.01)
canComp.SetLeftMargin(0.085)
canComp.SetTopMargin(0.11)
canComp.SetBottomMargin(0.11)

NfftPows = 6
NfftDers = 2
NfftBasis = NfftPows*NfftDers
#ufftBasisCRS = np.fromfile("../rotorBasis/fftRotorbasis.dat", dtype = np.double)
#ufftBasisCRS = np.resize(ufftBasisCRS, [NfftBasis, 40])
ufftBasis = np.fromfile("../movie/data/ufftRotorbasis.dat", dtype = np.double)
ufftBasis = np.resize(ufftBasis, [NfftBasis, 392])
lfftBasis = np.fromfile("../movie/data/lfftRotorbasis.dat", dtype = np.double)
lfftBasis = np.resize(lfftBasis, [NfftBasis, 604])
#for i in range(NfftBasis):
#  print(ufftBasis[i])
#print(ufftBasis[1])
#plt.plot(ufftBasis[0,:])
#plt.show()
#plt.plot(ufftBasis[1,:])
#plt.show()


tmp = np.fromfile("../movie/data/expValCos2_355.726000-360.026000.dat", dtype = np.double)
uCosEV = np.zeros((6,tmp.shape[0]), dtype = np.double)
uCosEV[0,:] = tmp
for i in range(1, uCosEV.shape[0]):
  uCosEV[i,:] = np.fromfile("../movie/data/expValCos" + str((i+1)*2) 
      + "_355.726000-360.026000.dat", dtype = np.double)
lCosEV = np.fromfile("/reg/d/psdm/amo/amoi0314/scratch/cosExpVals/expValCos2_357.045000-360.025000.dat", dtype = np.double)
print("lcosSize",lCosEV.shape)

Nlg = 5
Ncs = 5
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
      #print(ufftCoefRad[ilg,itr,:])
      itr += 1

NlRad = 30
lCosCoefRad = np.zeros((Nlg, Ncs, NlRad), dtype=np.double)
for ilg in range(Nlg):
  if (ilg>2):
    continue
  lInp = np.loadtxt("../movie/data/projections_172_e5_l" + str(2*ilg) 
     + ".dat", dtype = np.double)
  for i in range(Ncs):
    lCosCoefRad[ilg,i,:] = lInp[:,i]
    



####################
###  Make Plot  ###
####################

maxQ = 14.1207
Npix = 257
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


lCosCoef = np.zeros((Nlg, Ncs, Npix, Npix))
rInds = np.copy(dR)
rInds[rInds > Npix/2] = 0
rInds -= circRad*Npix/2.
rInds[rInds < 0] = 0
rInds *= 1.0*NlRad/((1 - circRad)*Npix/2.)
rIndsFlat_double = np.reshape(rInds, (-1, 1))
rIndsFlat = rIndsFlat_double.astype(np.int8)
for ilg in range(Nlg):
  for i in range(Ncs):
    lCosCoef[ilg,i,:,:] = np.reshape(lCosCoefRad[ilg,i,rIndsFlat], (Npix, Npix))

#plt.imshow(rInds)
#plt.show()


legImg = np.zeros((Nlg, Npix, Npix), dtype = np.double)
for ilg in range(Nlg):
  c = np.zeros(2*ilg + 1)
  c[-1] = 1
  legImg[ilg,:,:] = np.polynomial.legendre.legval(cosTheta, c)


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

ucosAvgs = [0.347042, 0.214022, 0.155776, 0.122814, 0.101518, 0.0865882]
lcosAvgs = [0.347042, 0.214022, 0.155776, 0.122814, 0.101518, 0.0865882]

uCSmax = np.argmax(ufftBasis[0])
uCSmin = np.argmin(ufftBasis[0][:uCSmax])
minV = 10
for i in range(uCSmin,uCSmax):
  if np.absolute(ufftBasis[0,i] - 0.3333333) < minV:
    minV = np.absolute(ufftBasis[0,i] - 0.3333333)
    uCSref = i 
#uCSref = np.argmin(np.absolute(ufftBasis[0][uCSmin:uCSmax] - 0.3333333))
lCSmax = np.argmax(lfftBasis[0])
lCSmin = np.argmin(lfftBasis[0,:lCSmax])
minV = 10
for i in range(lCSmin,lCSmax):
  if np.absolute(lfftBasis[0,i]) < minV:
    minV = np.absolute(lfftBasis[0,i])
    lCSref = i 

uCSref = 237
lCSref = 302
# Reference Image
uRefImg = np.zeros((Npix, Npix), dtype = np.double)
lRefImg = np.zeros((Npix, Npix), dtype = np.double)
for ilg in range(1,Nlg):
  for ibs in range(NfftBasis):
    uRefImg[:,:] += ufftCoef[ilg,ibs,:,:]*ufftBasis[ibs,uCSref]*legImg[ilg,:,:]
  for ics in range(Ncs):
    lRefImg[:,:] += lCosCoef[ilg,ics]*uCosEV[ics,uCSref]*legImg[ilg,:,:]


# Minimum Image
uMinImg = np.zeros((Npix, Npix), dtype = np.double)
lMinImg = np.zeros((Npix, Npix), dtype = np.double)
for ilg in range(1,Nlg):
  for ibs in range(NfftBasis):
    uMinImg[:,:] += ufftCoef[ilg,ibs,:,:]*ufftBasis[ibs,uCSmin]*legImg[ilg,:,:]
  for ics in range(Ncs):
    lMinImg[:,:] += lCosCoef[ilg,ics]*uCosEV[ics,uCSmin]*legImg[ilg,:,:]

# Maximum Image
uMaxImg = np.zeros((Npix, Npix), dtype = np.double)
lMaxImg = np.zeros((Npix, Npix), dtype = np.double)
for ilg in range(1,Nlg):
  for ibs in range(NfftBasis):
    uMaxImg[:,:] += ufftCoef[ilg,ibs,:,:]*ufftBasis[ibs,uCSmax]*legImg[ilg,:,:]
  for ics in range(Ncs):
    lMaxImg[:,:] += lCosCoef[ilg,ics]*uCosEV[ics,uCSmax]*legImg[ilg,:,:]

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

uTii = 355.726
uTff = 360.026
uTi = uTii + 2*(uTff - uTii)/44.
uTf = uTff - 2*(uTff - uTii)/44.
timeRatio = 1.75
timeWindow = 150
lOrigTi = 357.045
lOrigTf = 360.025
udt = (uTf - uTi)/ufftBasis[0].shape[0]
ldt = (lOrigTf - lOrigTi)/lfftBasis.shape[1]
uTref = udt*uCSref + 355.726
lTref = ldt*lCSref + 357.045
uCosHist = ROOT.TH1F("UEDcosSq", "UEDcosSq", \
    ufftBasis[0].shape[0], uTi, uTf)

lCosHist = ROOT.TH1F("LCLScosSq", "LCLScosSq", \
    lfftBasis.shape[1], 
    uTref - timeRatio*(uTref - (lOrigTi - (lTref - uTref))), 
    uTref + timeRatio*((lOrigTf - (lTref - uTref) - uTref)))

for i in range(ufftBasis[0].shape[0]):
  uCosHist.SetBinContent(i+1, ufftBasis[0,i])

for i in range(lfftBasis.shape[1]):
  lCosHist.SetBinContent(i+1, lfftBasis[0,i])

uCosHist.GetXaxis().SetRangeUser(uTref - timeWindow*udt, uTref + timeWindow*udt)


uColor = 4
lColor = 2
labelSize = 0.04
uCosHist.SetLineColor(uColor)
lCosHist.SetLineColor(lColor)
uCosHist.GetXaxis().SetTitleColor(uColor)
uCosHist.GetXaxis().SetLabelColor(uColor)
uCosHist.GetYaxis().SetTitle("#LT Cos^{2}(#theta) #GT")
uCosHist.GetXaxis().SetTitle("Time [ps]")
uCosHist.GetXaxis().CenterTitle(1)
uCosHist.GetYaxis().CenterTitle(1)

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

uCosHist.Draw()
canComp.Print("utest.png")
lCosHist.Draw()
canComp.Print("ltest.png")


uCosHist.Draw()
lCosHist.Draw("SAME")

cosLgnd = ROOT.TLegend(0.12, 0.7, 0.32, 0.82)
cosLgnd.AddEntry(uCosHist, "UED Simulation", "l")
cosLgnd.AddEntry(lCosHist, "LCLS Simulation", "l")
cosLgnd.Draw("SAME")

txt = ROOT.TText()
txt.SetTextSize(labelSize)
txt.SetTextAlign(22)
txt.SetTextColor(lColor)
axLabels = [357, 357.5, 358, 358.5, 359, 359.5]
for i in axLabels:
  txt.DrawText(i, uCosHist.GetMaximum()*1.18, str((i-uTref)/timeRatio + lTref)[:5])
txt.DrawText((uTref-timeWindow*udt + uTref+timeWindow*udt)/2, 
      uCosHist.GetMaximum()*1.33, "Time [ps]")


#plt.plot(l1dSimAln)
#plt.savefig("l1dSimAln")


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
txtLabel.SetTextSize(0.25)


pSize = 0.155
pSep = 0.004
p1x = 0.325
p1y = 0.46 
canComp.cd()
pad1a = ROOT.TPad("pad1a", "pad1a", p1x, p1y, p1x + pSize/2, p1y + pSize*xpix/ypix)
pad1a.SetRightMargin(0)
pad1a.SetLeftMargin(0)
pad1a.SetTopMargin(0)
pad1a.SetBottomMargin(0)
pad1a.SetFrameFillColor(0)
pad1a.SetFrameLineColor(0)
pad1a.SetBorderSize(0)
pad1a.Draw()
pad1a.cd()
lminHist.Draw("COL")
lcirc.Draw("SAME")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)
txtLabel.DrawText(0.0, 0.86, "A")

canComp.cd()
pad1b = ROOT.TPad("pad1b", "pad1b", p1x + pSize/2 + pSep, p1y, p1x + 2*pSize/2 + pSep, p1y + pSize*xpix/ypix)
pad1b.SetRightMargin(0)
pad1b.SetLeftMargin(0)
pad1b.SetTopMargin(0)
pad1b.SetBottomMargin(0)
pad1b.SetFrameFillColor(0)
pad1b.SetFrameLineColor(0)
pad1b.SetBorderSize(0)
pad1b.Draw()
pad1b.cd()
uminHist.Draw("COL")
ucirc.Draw("SAME")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)


p2x = 0.56
p2y = 0.15
canComp.cd()
pad2a = ROOT.TPad("pad2a", "pad2a", p2x, p2y, p2x + pSize/2, p2y + pSize*xpix/ypix)
pad2a.SetRightMargin(0)
pad2a.SetLeftMargin(0)
pad2a.SetTopMargin(0)
pad2a.SetBottomMargin(0)
pad2a.SetFrameFillColor(0)
pad2a.SetBorderSize(0)
pad2a.Draw()
pad2a.cd()
lmaxHist.Draw("COL")
txt.SetTextSize(labelSize*1.2)
lcirc.Draw("SAME")
txtLabel.DrawText(0.0, 0.86, "B")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)

canComp.cd()
pad2b = ROOT.TPad("pad2b", "pad2b", p2x + pSize/2 + pSep, p2y, p2x + 2*pSize/2 + pSep, p2y + pSize*xpix/ypix)
pad2b.SetRightMargin(0)
pad2b.SetLeftMargin(0)
pad2b.SetTopMargin(0)
pad2b.SetBottomMargin(0)
pad2b.SetFrameFillColor(0)
pad2b.SetBorderSize(0)
pad2b.Draw()
pad2b.cd()
umaxHist.Draw("COL")
ucirc.Draw("SAME")
lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)


canComp.cd()
canComp.Update()
canComp.Print("compareTimeBasis.png")




#####################
##  Compare Coeff  ##
#####################

uCoefImg = np.zeros((Npix, Npix), dtype = np.double)
lCoefImg = np.zeros((Npix, Npix), dtype = np.double)
for ilg in range(1,Nlg):
  uCoefImg[:,:] += ufftCoef[ilg,0,:,:]*legImg[ilg,:,:]
  lCoefImg[:,:] += lCosCoef[ilg,0]*legImg[ilg,:,:]

uNorm = 1.0/max(np.absolute([np.amax(uCoefImg), np.amin(uCoefImg)]))
lNorm = 1.0/max(np.absolute([np.amax(lCoefImg), np.amin(lCoefImg)]))

uCoefImg *= uNorm
lCoefImg *= lNorm




######################
##  Power Spectrum  ##
######################


uPowSpect = np.zeros((NfftBasis,1), dtype=np.double)
lPowSpect = np.zeros((Ncs,1), dtype=np.double)
for ilg in range(Nlg):
  for i in range(NuRad):
    for ibs in range(NfftBasis):
      uPowSpect[ibs] += ufftCoefRad[ilg, ibs, i]**2
    for ics in range(Ncs):
      if ilg < 2 and i < NlRad:
        lPowSpect[ics] += lCosCoefRad[ilg, ics, i]**2
print(uPowSpect)
print(sum(uPowSpect))
uPowSpect /= sum(uPowSpect)
lPowSpect /= sum(lPowSpect)

uPowHist = ROOT.TH1F("uPowHist", "uPowHist", Ncs, 1, Ncs+1)
lPowHist = ROOT.TH1F("lPowHist", "pPowHist", Ncs, 1, Ncs+1)
for i in range(Ncs):
  uPowHist.SetBinContent(i+1, uPowSpect[i])
  lPowHist.SetBinContent(i+1, lPowSpect[i])
  uPowHist.SetBinError(i+1, 1e-10)
  lPowHist.SetBinError(i+1, 1e-10)

uPowHist.GetXaxis().SetTitle("Time Basis")
uPowHist.GetYaxis().SetTitle("Power")
uPowHist.GetYaxis().SetTitleSize(0.13)
uPowHist.GetXaxis().SetTitleSize(0.13)
uPowHist.GetXaxis().SetLabelSize(0)
uPowHist.GetYaxis().SetLabelSize(labelSize*1.7)
uPowHist.GetYaxis().SetTitleOffset(0.33)
uPowHist.GetXaxis().SetTitleOffset(0.77)
uPowHist.GetYaxis().CenterTitle(1)
uPowHist.GetXaxis().CenterTitle(1)
uPowHist.SetMarkerColor(uColor)
uPowHist.SetMarkerSize(1.6)
uPowHist.SetLineColor(uColor)
lPowHist.SetMarkerColor(lColor)
lPowHist.SetMarkerSize(1.6)
lPowHist.SetLineColor(lColor)

powLgnd = ROOT.TLegend(0.67, 0.78, 0.875, 0.95)
powLgnd.SetBorderSize(1)
powLgnd.AddEntry(uPowHist, "UED Simulation", "lp")
powLgnd.AddEntry(lPowHist, "LCLS Simulation", "lp")

canCoeff.cd()
ccpad3 = ROOT.TPad("ccpad3", "ccpad3", 0, 0, 1, (ycut-20)/ypix)
ccpad3.SetRightMargin(0.09)
ccpad3.SetLeftMargin(0.09)
ccpad3.SetTopMargin(0)
ccpad3.SetBottomMargin(0.25)
ccpad3.SetFrameFillColor(0)
ccpad3.SetFrameLineColor(0)
ccpad3.SetBorderSize(0)
ccpad3.SetLogy(1)
ccpad3.Draw()
ccpad3.cd()


uPowHist.Draw("LPE")
lPowHist.Draw("SAMELPE")
powLgnd.Draw("SAME")

tex = ROOT.TLatex()
tex.SetTextSize(0.08)
axLabels = ["Basis 1", "Basis 2", "Basis 3", "Basis 4", "Basis 5"]

print(uPowHist.GetMinimum())
for i in range(len(axLabels)):
  #tex.DrawLatex(1.25 + i, -1, axLabels[i])
  tex.DrawLatexNDC(0.14 + 0.165*i, 0.16, axLabels[i])

canCoeff.Print("compareCoeff.png")

"""

