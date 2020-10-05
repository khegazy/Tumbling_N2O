import sys, glob
sys.argv.append("-b")
import ROOT 
import numpy as np
import pickle as pl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = u'/reg/neh/home5/khegazy/analysis/tumblingN2O/movie/env/bin/ffmpeg'

#import root_numpy

#from root_numpy import hist2array


sys.path.insert(0, "/reg/neh/home/khegazy/baseScripts/")
from pyPlotFunctions import PLOTclass 
#from plotGstyle import genStyle, gen2dStyle


plc = PLOTclass()
ROOT.gROOT.SetStyle("genStyle")
xpix = 1500
ypix = 600
ypixRat = 600.0/ypix
canComp = ROOT.TCanvas("canComp", "canComp", xpix, ypix);

rMarg = 0.01
lMarg = 0.085
tMarg = 0.01
bMarg = 0.11
canComp.SetRightMargin(rMarg)
canComp.SetLeftMargin(lMarg)
canComp.SetTopMargin(tMarg)
canComp.SetBottomMargin(bMarg)

NfftPows = 6
NfftDers = 2
NfftBasis = NfftPows*NfftDers
#ufftBasisCRS = np.fromfile("../rotorBasis/fftRotorbasis.dat", dtype = np.double)
#ufftBasisCRS = np.resize(ufftBasisCRS, [NfftBasis, 40])
ufftBasis = np.fromfile("../rotorBasis/uedfftRotorbasis.dat", dtype = np.double)
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


#uCosEV = np.fromfile("../movie/data/expValCos2_355.726000-360.026000.dat", dtype = np.double)
uCosEV = np.fromfile("../simulation/rotation/output/expVals/expValCos2_315.914000-320.214000-0.010000.dat", dtype = np.double)
lCosEV = np.fromfile("/reg/d/psdm/amo/amoi0314/scratch/cosExpVals/expValCos2_317.000000-320.500000-0.004230.dat", dtype = np.double)
uIndMax = np.argmax(uCosEV)
uIndMin = np.argmin(uCosEV)
uIndRef = uIndMin + np.argmin(np.abs(uCosEV[uIndMin:uIndMax]-0.333333))
evMax = np.argmax(lCosEV)
evMin = np.argmin(lCosEV[:evMax])
lIndRef = evMin + np.argmin(np.abs(lCosEV[evMin:evMax]-0.333333))

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
# Sum bin 10 - 14
print("Getting proj coef")
lfftCoefRad = np.zeros((Nlg, NfftBasis, NlRad), dtype=np.double)
for ilg in range(Nlg):
  #if (ilg>2):
  #  continue
  lInp = np.loadtxt("../movie/data/projections-amoi0314-r0172_e5_g1_l" + str(2*ilg) 
     + "_tspan1.083.dat", dtype = np.double)
  for i in range(NfftBasis):
    lfftCoefRad[ilg,i,:] = lInp[:,i+1]
    
print("got proj coef")



####################
###  Make Plot  ###
####################

fig = plt.figure()
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
ax = fig.add_subplot(111)
for ilg in range(Nlg):
  for i in range(NfftBasis):
    lfftCoef[ilg,i,:,:] = np.reshape(lfftCoefRad[ilg,i,rIndsFlat], (Npix, Npix))
    ax.clear()
    ax.plot(lfftCoefRad[ilg,i,:])
    plt.savefig("lfftBasisCoeff_basis-"+str(i)+"_lg-"+str(ilg)+".png")

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

plotAngDist = False 
if plotAngDist:
  ind = 0
  for f in glob.glob("../simulation/rotation/output/alignmentPDFs/UED/plotJobs/*dat"):
  #for f in glob.glob("../simulation/rotation/output/alignmentPDFs/UED/test/*dat"):
    print("file",f)

    pdf = np.fromfile(f, dtype=np.double)
    thBins = 200
    phBins = 200
    print("orig pdf shape",pdf.shape)
    pdf = np.reshape(pdf, (thBins,phBins))
    #plt.imshow(pdf)
    #plt.show()
    deltaTh = 3.14159265359/thBins
    deltaPh = 2*3.14159265359/phBins
    grph = ROOT.TGraph2D()
    ind = 0

    #plt.imshow(pdf)
    #plt.show()
    #continue

    ax = fig.add_subplot(111, projection="3d")
    Xg, Yg = np.mgrid[0:thBins+1,0:phBins+1]
    Xg[-1,:] = Xg[0,:]
    Yg[:,-1] = Yg[:,0]
    r = 1e6*pdf[Xg, Yg]
    theta = Xg*deltaTh 
    phi = Yg*deltaPh + 3.14159
    X = r*np.sin(theta)*np.sin(phi)
    Y = r*np.sin(theta)*np.cos(phi)
    #print("X")
    #print(Xg)
    #print("rows")
    #print(pdf[:,50])
    #print("cols")
    #print(pdf[50,:])
    Z = r*np.cos(theta)
    Xf = np.zeros((2*Xg.shape[0]-1, Xg.shape[1]))
    Yf = np.zeros((2*Xg.shape[0]-1, Xg.shape[1]))
    Zf = np.zeros((2*Xg.shape[0]-1, Xg.shape[1]))
    #print(np.sum(r))
    #exit(0)
    for i in range(X.shape[0]-1):
      Xf[i,:] = X[i,:]
      Xf[-i-1,:] = X[i,:]
      Yf[i,:] = Y[i,:]
      Yf[-i-1,:] = Y[i,:]
      Zf[i,:] = Z[i,:]
      Zf[-i-1,:] = -1*Z[i,:]
    surf = ax.plot_surface(Xf,Yf,Zf, rstride=8, cstride=8, cmap="jet", linewidth=0.4)
    #surf = ax.plot_surface(Xf,Yf,Zf, rstride=2, cstride=2, cmap="jet")
    ax.view_init(azim=0,elev=0)
    #Z = -1*Z
    #surf = ax.plot_surface(X,Y,Z, rstride=5, cstride=5)

    print("drawing")
    #grph.Draw("SURF1")
    #canComp.Print("testingB.png")

    #surf = ax.plot_surface(X,Y,Z, color='b')
    #surf = ax.scatter(X,Y,Z, color='b')
    #surf = ax.plot_trisurf(X,Y,Z, linewidth=0.2, antialiased=True)
    ax.set_xlim([-50,50])
    ax.set_ylim([-50,50])
    ax.set_zlim([-50,50])
    plt.axis('off')
    #fig.colorbar(surf)
    #plt.show()
    plt.savefig(f[f.find("Time"):-4]+".png")
    print("drawn")





########################
uCSref = 237
lCSref = 302
uTi = 316 #355.726
uTf = 320.5 #360.026
timeRatio = 1.75
indTimeWindow = 150
lOrigTi = 357.045
lOrigTf = 360.025
udt = (uTf - uTi)/ufftBasis[0].shape[0]
ldt = (lOrigTf - lOrigTi)/lfftBasis.shape[1]
uTmax = udt*uIndMax + uTi
uTmin = udt*uIndMin + uTi
uTref = udt*uIndRef + uTi
lTref = ldt*lIndRef + 357.045
uCosHist = ROOT.TH1F("UEDcosSq", "UEDcosSq", \
    ufftBasis[0].shape[0], uTi, uTf)
lShift = 2*timeRatio*(uTref - (lOrigTi - (lTref - uTref)))*(lCSref - lfftBasis.shape[1]/2.)/float(lfftBasis.shape[1])
lCosHist = ROOT.TH1F("LCLScosSq", "LCLScosSq", \
    lfftBasis.shape[1], 
    lShift + uTref - timeRatio*(uTref - (lOrigTi - (lTref - uTref))), 
    lShift + uTref + timeRatio*((lOrigTf - (lTref - uTref) - uTref)))

for i in range(ufftBasis[0].shape[0]):
  uCosHist.SetBinContent(i+1, uCosEV[i])

for i in range(lfftBasis.shape[1]):
  lCosHist.SetBinContent(i+1, lCosEV[i])

tMin = uTref - indTimeWindow*udt
tMax = uTref + indTimeWindow*udt
timeWindow = tMax - tMin
print("limits", timeWindow*udt, uTref - timeWindow*udt, uTref + timeWindow*udt)
uCosHist.GetXaxis().SetRangeUser(tMin, tMax)
yMin = 0.22 
uCosHist.SetMinimum(yMin)
yMax = 0.47 
uCosHist.SetMaximum(yMax)
yWindow = yMax - yMin

uColor = 4
lColor = 2
labelSize = 0.04
uCosHist.SetLineColor(uColor)
uCosHist.SetLineWidth(3)
lCosHist.SetLineColor(lColor)
#uCosHist.GetXaxis().SetTitleColor(uColor)
#uCosHist.GetXaxis().SetLabelColor(uColor)
uCosHist.GetYaxis().SetTitle("#LT cos^{2}(#theta) #GT")
uCosHist.GetXaxis().SetTitle("Time [ps]")
uCosHist.GetXaxis().CenterTitle(1)
uCosHist.GetYaxis().CenterTitle(1)

uCosHist.Draw("C")

txt = ROOT.TText()
txt.SetTextSize(labelSize)
txt.SetTextAlign(22)
txt.SetTextColor(lColor)
axLabels = [357, 357.5, 358, 358.5, 359, 359.5]
#for i in axLabels:
#  txt.DrawText(i, uCosHist.GetMaximum()*1.18, str((i-uTref)/timeRatio + lTref)[:5])
#txt.DrawText((uTref-timeWindow*udt + uTref+timeWindow*udt)/2, 
#      uCosHist.GetMaximum()*1.33, "Time [ps]")



#plt.plot(ufftBasis[0,:])
#plt.plot(ufftBasis[2,:])
#plt.plot(ufftBasis[4,:])
#plt.show()

test = np.zeros(5)
print("test",test)

Ninsets = 5
ltime = []
deltaT = (float(uTmax) - float(uTmin))/3
for i in range(30):
  tm = float(uTmin) - 10*deltaT + i*deltaT
  if (tm>tMin) and (tm<tMax):
    if abs(tm-uTmax)>0.1 and abs(tm-uTmin)>0.1:
      ltime.append(tm)

time = np.array(ltime)
px = (time - tMin)/timeWindow
px[-1] -=0.01

print("tmax and tmin", uTmax, uTmin)
print("evaluation times", time)
py = bMarg + (1-bMarg-tMarg)*(uCosEV[(((tMin + px*timeWindow)-uTi)/udt).astype(int)] - yMin)/yWindow
print("py", py)
py[py>1.0-0.155*xpix/ypix] = 0.95 - 0.155*xpix/ypix

lImgs = np.zeros((len(px),Npix, Npix), dtype = np.double)
uImgs = np.zeros((len(px),Npix, Npix), dtype = np.double)
for i in range(len(px)):
  timeInd = int((time[i] - uTi)/udt)

  for ilg in range(Nlg):
    for ibs in range(NfftBasis):
      lImgs[i,:,:] += lfftCoef[ilg,ibs]*ufftBasis[ibs,timeInd]*legImg[ilg,:,:]
      if ilg is not 0:
        uImgs[i,:,:] += ufftCoef[ilg,ibs]*ufftBasis[ibs,timeInd]*legImg[ilg,:,:]

  # NEED THIS or else first two images are blank?????
  if i<2 and True:
    plt.imshow(lImgs[i,:,:])
    #plt.colorbar()
    #plt.show()
    plt.imshow(uImgs[i,:,:])
    #plt.colorbar()
    #plt.show()

lMaxImg = np.zeros((Npix, Npix), dtype = np.double)
lMinImg = np.zeros((Npix, Npix), dtype = np.double)
uMaxImg = np.zeros((Npix, Npix), dtype = np.double)
uMinImg = np.zeros((Npix, Npix), dtype = np.double)
for ilg in range(Nlg):
  for ibs in range(NfftBasis):
    lMinImg[:,:] += lfftCoef[ilg,ibs]*ufftBasis[ibs,int((uTmin - uTi)/udt)]*legImg[ilg,:,:]
    lMaxImg[:,:] += lfftCoef[ilg,ibs]*ufftBasis[ibs,int((uTmax - uTi)/udt)]*legImg[ilg,:,:]
    if ilg is not 0:
      uMinImg[:,:] += ufftCoef[ilg,ibs]*ufftBasis[ibs,int((uTmin - uTi)/udt)]*legImg[ilg,:,:]
      uMaxImg[:,:] += ufftCoef[ilg,ibs]*ufftBasis[ibs,int((uTmax - uTi)/udt)]*legImg[ilg,:,:]


lNorm = 1.0/max(np.absolute([np.amax(lMaxImg), np.amin(lMaxImg), np.amax(lMinImg), np.amin(lMinImg), np.amax(lImgs), np.amin(lImgs)]))
uNorm = 1.0/max(np.absolute([np.amax(uMaxImg), np.amin(uMaxImg), np.amax(uMinImg), np.amin(uMinImg), np.amax(uImgs), np.amin(uImgs)]))
#lNorm = 1.0/max(np.absolute([np.amax(lImgs)]))
lImgs *= lNorm
lMaxImg *= lNorm
lMinImg *= lNorm
uImgs *= uNorm
uMaxImg *= uNorm
uMinImg *= uNorm
print("maxmax", np.amax(lMaxImg))

line = ROOT.TLine()
line.SetLineWidth(1)
lb = ROOT.TLine()
rb = ROOT.TLine()
tb = ROOT.TLine()
bb = ROOT.TLine()
lb.SetLineColor(0)
rb.SetLineColor(0)
tb.SetLineColor(0)
bb.SetLineColor(0)

circRad = 0.05
lcirc = ROOT.TEllipse(1.0, 0.5, 2*circRad, circRad, 90, 270)
lcirc.SetFillColor(1)
ucirc = ROOT.TEllipse(0.0, 0.5, 2*circRad, circRad)
ucirc.SetFillColor(1)


pSize = 0.155
pSep = 0.001
###############################
#####  Max and Min Hists  #####
###############################

lminHist = ROOT.TH2F("lminHist", "lminHist", Npix/2+1, 0, 1, Npix, 0, 1)
lmaxHist = ROOT.TH2F("lmaxHist", "lmaxHist", Npix/2+1, 0, 1, Npix, 0, 1)
uminHist = ROOT.TH2F("uminHist", "uminHist", Npix/2+1, 0, 1, Npix, 0, 1)
umaxHist = ROOT.TH2F("umaxHist", "umaxHist", Npix/2+1, 0, 1, Npix, 0, 1)
lminHist.SetContour(80)
lmaxHist.SetContour(80)
uminHist.SetContour(80)
umaxHist.SetContour(80)
for ir in range(Npix):
  for ic in range(Npix/2 + 1):
    lminHist.SetBinContent(ic+1, ir+1, lMinImg[ir,ic])
    lmaxHist.SetBinContent(ic+1, ir+1, lMaxImg[ir,ic])
    uminHist.SetBinContent(ic+1, ir+1, uMinImg[ir,Npix/2 + ic])
    umaxHist.SetBinContent(ic+1, ir+1, uMaxImg[ir,Npix/2 + ic])

for ir in range(Npix):
  for ic in range(Npix/2 + 1):
    if np.sqrt((ir-Npix/2)**2 + (ic-Npix/2)**2) > Npix/2:
      lmaxHist.SetBinContent(ic+1, ir+1, -1e20)
      lminHist.SetBinContent(ic+1, ir+1, -1e20)
      umaxHist.SetBinContent(Npix/2 - ic + 2, ir+1, -1e20)
      uminHist.SetBinContent(Npix/2 - ic + 2, ir+1, -1e20)


lmaxHist.SetMaximum(1)
lminHist.SetMaximum(1)
umaxHist.SetMaximum(1)
uminHist.SetMaximum(1)
lmaxHist.SetMinimum(-1)
lminHist.SetMinimum(-1)
umaxHist.SetMinimum(-1)
uminHist.SetMinimum(-1)
lmaxHist.GetXaxis().SetLabelSize(0)
lminHist.GetXaxis().SetLabelSize(0)
umaxHist.GetXaxis().SetLabelSize(0)
uminHist.GetXaxis().SetLabelSize(0)
lmaxHist.GetXaxis().SetLabelOffset(999)
lminHist.GetXaxis().SetLabelOffset(999)
umaxHist.GetXaxis().SetLabelOffset(999)
uminHist.GetXaxis().SetLabelOffset(999)
lmaxHist.GetYaxis().SetLabelSize(0)
lminHist.GetYaxis().SetLabelSize(0)
umaxHist.GetYaxis().SetLabelSize(0)
uminHist.GetYaxis().SetLabelSize(0)
lmaxHist.GetXaxis().SetTickSize(0)
lminHist.GetXaxis().SetTickSize(0)
umaxHist.GetXaxis().SetTickSize(0)
uminHist.GetXaxis().SetTickSize(0)
lmaxHist.GetYaxis().SetTickSize(0)
lminHist.GetYaxis().SetTickSize(0)
umaxHist.GetYaxis().SetTickSize(0)
uminHist.GetYaxis().SetTickSize(0)


#pxM = 0.6
#pyM = 0.16
pxM = 0.585
pyM = 0.38
canComp.cd()
padM = ROOT.TPad("padM", "padM", pxM, pyM, pxM + pSize/2, pyM + pSize*xpix/ypix)
#pad = ROOT.TPad("pad", "pad", p1x, p1y, p1x + pSize/2, p1y + pSize*xpix/ypix)
padM.SetRightMargin(0)
padM.SetLeftMargin(0)
padM.SetTopMargin(0)
padM.SetBottomMargin(0)
padM.SetFrameFillColor(0)
padM.SetFrameLineColor(0)
padM.SetBorderSize(0)
padM.Draw()
padM.cd()

lmaxHist.Draw("COL")
lcirc.Draw("SAME")

lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)
#txtLabel.DrawText(0.0, 0.86, "A")

canComp.cd()
line.DrawLineNDC(pxM + pSize/2 + pSep, pyM + pSize*xpix/ypix,\
    lMarg + (1-rMarg-lMarg)*(uTmax - tMin)/timeWindow,\
    bMarg + (1-bMarg-tMarg)*(uCosEV[int((uTmax - uTi)/udt)] - yMin)/yWindow)

canComp.cd()
padM1 = ROOT.TPad("padM1", "padM1", pSep*2 + pxM + pSize/2, pyM, pSep*2 + pxM + pSize, pyM + pSize*xpix/ypix)
padM1.SetRightMargin(0)
padM1.SetLeftMargin(0)
padM1.SetTopMargin(0)
padM1.SetBottomMargin(0)
padM1.SetFrameFillColor(0)
padM1.SetFrameLineColor(0)
padM1.SetBorderSize(0)
padM1.Draw()
padM1.cd()

umaxHist.Draw("COL")
ucirc.Draw("SAME")

lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)



#pxm = 0.35
#pym = 0.6
pxm = 0.345
pym = 0.28
canComp.cd()
padm = ROOT.TPad("padm", "padm", pxm, pym, pxm + pSize/2, pym + pSize*xpix/ypix)
#pad = ROOT.TPad("pad", "pad", p1x, p1y, p1x + pSize/2, p1y + pSize*xpix/ypix)
padm.SetRightMargin(0)
padm.SetLeftMargin(0)
padm.SetTopMargin(0)
padm.SetBottomMargin(0)
padm.SetFrameFillColor(0)
padm.SetFrameLineColor(0)
padm.SetBorderSize(0)
padm.Draw()
padm.cd()

lminHist.Draw("COL")
lcirc.Draw("SAME")

lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)
#txtLabel.DrawText(0.0, 0.86, "A")

canComp.cd()
line.DrawLineNDC(pxm + pSize/2 + pSep, pym,\
    lMarg + (1-rMarg-lMarg)*(uTmin - tMin)/timeWindow,\
    bMarg + (1-bMarg-tMarg)*(uCosEV[int((uTmin - uTi)/udt)] - yMin)/yWindow)

canComp.cd()
padm1 = ROOT.TPad("padm1", "padm1", pSep*2 + pxm + pSize/2, pym, pSep*2 + pxm + pSize, pym + pSize*xpix/ypix)
padm1.SetRightMargin(0)
padm1.SetLeftMargin(0)
padm1.SetTopMargin(0)
padm1.SetBottomMargin(0)
padm1.SetFrameFillColor(0)
padm1.SetFrameLineColor(0)
padm1.SetBorderSize(0)
padm1.Draw()
padm1.cd()

uminHist.Draw("COL")
ucirc.Draw("SAME")

lb.DrawLineNDC(0, 0, 0, 1)
rb.DrawLineNDC(1, 0, 1, 1)
tb.DrawLineNDC(0, 1, 1, 1)
bb.DrawLineNDC(0, 0, 1, 0)

#############################
#####  All other hists  #####
#############################

pSize = 0.155/2.
ofs = 0.17
alt = np.arange(Ninsets, dtype=float)
alt[alt%2==1] = -1
alt[alt>=0] = 1
yOff = alt*ofs
print(yOff)
#py -= 0.08 
#py[py+pSize*xpix/ypix>0.98] -= 2*ofs
#py[py<0.06] += 2*ofs

px[0] += 0.07
py[0] += 0.12
px[1] += 0.01
py[1] -= 0.3
px[2] += 0.03
py[2] += 0.09
px[3] -= 0.05
py[3] -= 0.16
px[4] += 0.045
py[4] -= 0.17
px[5] -= 0.04
py[5] += 0.17
px[6] += 0.035
py[6] += 0.18
px[7] -= 0.03
py[7] -= 0.27
px[8] -= 0.045
py[8] += 0.1
px[9] -= 0.07
py[9] -= 0.27

print("change py")
#py = np.ones((Ninsets))*0.5
lnSX = px + pSize/2
lnSY = py + pSize*xpix/ypix
lnEX = lMarg + (1-rMarg-lMarg)*(time - tMin)/timeWindow
lnEY = bMarg + (1-bMarg-tMarg)*(uCosEV[((time - uTi)/udt).astype(int)] - yMin)/yWindow
lnSY[0] = py[0]
lnSY[2] = py[2]
lnSX[3] = px[3] + pSep + pSize
lnSY[3] = py[3] + pSize*xpix/ypix/2
lnSX[4] = px[4]
lnSY[4] = py[4] + pSize*xpix/ypix/2
lnSY[5] = py[5]
lnSX[6] = px[6]
lnSY[6] = py[6] + pSize*xpix/ypix/2
lnSY[8] = py[8]
#lnSY[10] = py[10]
#lnSY[12] = py[12]
#lnSY[14] = py[14]
#line.DrawLineNDC(pxM + pSize/2, pyM + pSize*xpix/ypix,\
#    lMarg + (1-rMarg-lMarg)*(uTmax - tMin)/timeWindow,\
#    bMarg + (1-bMarg-tMarg)*(ufftBasis[0,int((uTmax - uTi)/udt)] - yMin)/yWindow)
lHists = []
uHists = []
#lHist = ROOT.TH2F("lminHist", "lminHist", Npix/2+1, 0, 1, Npix, 0, 1)
#lHist.SetContour(80)
print("start loop")
for i in range(len(px)):
  lHists.append(ROOT.TH2F("lHist"+str(i), "lHist"+str(i), Npix/2+1, 0, 1, Npix, 0, 1))
  uHists.append(ROOT.TH2F("uHist"+str(i), "uHist"+str(i), Npix/2+1, 0, 1, Npix, 0, 1))
  lHists[i].SetContour(80)
  uHists[i].SetContour(80)
  for ir in range(Npix):
    for ic in range(Npix/2 + 1):
      lHists[i].SetBinContent(ic+1, ir+1, lImgs[i,ir,ic])
      uHists[i].SetBinContent(ic+1, ir+1, uImgs[i,ir,Npix/2 + ic])

  for ir in range(Npix):
    for ic in range(Npix/2 + 1):
      if np.sqrt((ir-Npix/2)**2 + (ic-Npix/2)**2) > Npix/2:
        lHists[i].SetBinContent(ic+1, ir+1, -1e20)
        uHists[i].SetBinContent(Npix/2 - ic + 2, ir+1, -1e20)

  print(lHists[i].GetMaximum(), lHists[i].GetMinimum())
  lHists[i].SetMaximum(1)
  uHists[i].SetMaximum(1)
  lHists[i].SetMinimum(-1)
  uHists[i].SetMinimum(-1)
  lHists[i].GetXaxis().SetLabelSize(0)
  uHists[i].GetXaxis().SetLabelSize(0)
  lHists[i].GetXaxis().SetLabelOffset(999)
  uHists[i].GetXaxis().SetLabelOffset(999)
  lHists[i].GetYaxis().SetLabelSize(0)
  uHists[i].GetYaxis().SetLabelSize(0)
  lHists[i].GetXaxis().SetTickSize(0)
  uHists[i].GetXaxis().SetTickSize(0)
  lHists[i].GetYaxis().SetTickSize(0)
  uHists[i].GetYaxis().SetTickSize(0)


  canComp.cd()
  pad = ROOT.TPad("pad", "pad", px[i], py[i], px[i] + pSize/2, py[i] + pSize*xpix/ypix)
  #pad = ROOT.TPad("pad", "pad", p1x, p1y, p1x + pSize/2, p1y + pSize*xpix/ypix)
  pad.SetRightMargin(0)
  pad.SetLeftMargin(0)
  pad.SetTopMargin(0)
  pad.SetBottomMargin(0)
  pad.SetFrameFillColor(0)
  pad.SetFrameLineColor(0)
  pad.SetBorderSize(0)
  pad.Draw()
  pad.cd()
  
  lHists[i].Draw("COL")
  lcirc.Draw("SAME")

  lb.DrawLineNDC(0, 0, 0, 1)
  rb.DrawLineNDC(1, 0, 1, 1)
  tb.DrawLineNDC(0, 1, 1, 1)
  bb.DrawLineNDC(0, 0, 1, 0)
  #txtLabel.DrawText(0.0, 0.86, "A")

  canComp.cd()
  pad1 = ROOT.TPad("pad1", "pad1", px[i] + pSep + pSize/2, py[i], px[i] + pSep + pSize, py[i] + pSize*xpix/ypix)
  pad1.SetRightMargin(0)
  pad1.SetLeftMargin(0)
  pad1.SetTopMargin(0)
  pad1.SetBottomMargin(0)
  pad1.SetFrameFillColor(0)
  pad1.SetFrameLineColor(0)
  pad1.SetBorderSize(0)
  pad1.Draw()
  pad1.cd()
  
  uHists[i].Draw("COL")
  ucirc.Draw("SAME")

  lb.DrawLineNDC(0, 0, 0, 1)
  rb.DrawLineNDC(1, 0, 1, 1)
  tb.DrawLineNDC(0, 1, 1, 1)
  bb.DrawLineNDC(0, 0, 1, 0)

  canComp.cd()
  line.DrawLineNDC(lnSX[i] + pSep/2., lnSY[i], lnEX[i], lnEY[i])

print("finished loop")


canComp.cd()
uCosHist.Draw("CSAME")
canComp.Update()
canComp.Print("uedFrameData.png")



"""
pnc bank
channel services
500 first ave
pittsburg pa, 15219
mailstop p7-pfsc-01-s
"""

