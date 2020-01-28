import sys
sys.argv.append("-b")
import ROOT 
import numpy as np
import matplotlib.pyplot as plt
#import root_numpy
#from root_numpy import hist2array

sys.path.insert(0, "/reg/neh/home/khegazy/baseTools/tools/")
from pyPlotFunctions import PLOTclass 
#from plotGstyle import genStyle, gen2dStyle


plc = PLOTclass()
ROOT.gROOT.SetStyle("genStyle")
ROOT.gStyle.SetErrorX(0)

mainDir = "/reg/neh/home/khegazy/analysis/tumblingN2O/UED/output/"
uDataFile = ROOT.TFile.Open(mainDir + "data/combinedResults.root", "read")
uSimRegFile = ROOT.TFile.Open(mainDir + "MC/bend2/combinedResults.root", "read")

lData = np.zeros((30, 206), dtype=float)
inpLfileNames = ["/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/fitParameters/data/diff_legendre-amoi0314-r0172_e3_g1_l2_tspan1.263.dat","/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/fitParameters/data/diff_legendre-amoi0314-r0172_e4_g1_l2_tspan1.263.dat", "/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/fitParameters/data/diff_legendre-amoi0314-r0172_e5_g1_l2_tspan1.263.dat","/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/fitParameters/data/diff_legendre-amoi0314-r0172_e6_g1_l2_tspan1.263.dat", "/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/fitParameters/data/diff_legendre-amoi0314-r0172_e7_g1_l2_tspan1.263.dat"];
for flName in inpLfileNames :
  fl = open(flName, 'r')
  data = np.loadtxt(fl, dtype=np.str, delimiter='\t')
  np.delete(data, 0, 0)
  data = data[:,:-1]
  data = data[:,14:-36]
  data = data[:,:].astype(np.float)
  lData += data

mainSimDir = "/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/fitParameters/output/"
u1dData = np.fromfile(mainSimDir + "UEDdata.dat", dtype=np.double)
u1dDataStd = np.sqrt(np.fromfile(mainSimDir + "UEDdatVar.dat", dtype=np.double))
u1dSclData = np.fromfile(mainSimDir + "UEDscldata.dat", dtype=np.double)
u1dSclSim = np.fromfile(mainSimDir + "UEDsclSim.dat", dtype=np.double)
u1dSclSim = u1dSclSim[:u1dData.size]
u1dSimAln = np.fromfile(mainSimDir + "UEDsimAln.dat", dtype=np.double)
l1dData = np.fromfile(mainSimDir + "LCLSdata.dat", dtype=np.double)
l1dDataStd = np.sqrt(np.fromfile(mainSimDir + "LCLSdataVar.dat", dtype=np.double))
l1dSclData = np.fromfile(mainSimDir + "LCLSsclData.dat", dtype=np.double)
l1dSim = np.fromfile(mainSimDir + "LCLSsim.dat", dtype=np.double)
l1dSclSim = np.fromfile(mainSimDir + "LCLSsclSim.dat", dtype=np.double)
l1dSimAln = np.fromfile(mainSimDir + "LCLSsimAln.dat", dtype=np.double)

txt = ROOT.TText()
txt.SetTextSize(0.08)
txt.SetTextAlign(22)
txt.SetTextFont(62)



#plt.plot(u1dSimAln)
#plt.savefig("u1dSimAln")
#plt.plot(l1dSimAln)
#plt.savefig("l1dSimAln")

#####  Compare legendres and fit  #####

maxRange = 1 

dNorm = 0
sNorm = 0

uDataHistInp = uDataFile.Get("legendre2Comb").Clone("udHclone")
uDataHist = ROOT.TH2F("UEDdata", "UEDdata", uDataHistInp.GetNbinsX(), 0, 1, uDataHistInp.GetNbinsY(), 0, 14.1207)
uDataHist.GetXaxis().SetLabelSize(0)
lDataHist = ROOT.TH2F("LCLSdata", "LCLSdata", lData.shape[1], 1.7125, 2.9796, lData.shape[0], 500, 530)
lDataHist.GetXaxis().SetLabelSize(0)
lDataHist.SetContour(80)


for ix in range(1, uDataHistInp.GetNbinsX() + 1):
  for iy in range(1, uDataHist.GetNbinsY() + 1):
    uDataHist.SetBinContent(ix, iy, uDataHistInp.GetBinContent(ix, iy))

for ir in range(lData.shape[0]) :
  for ic in range(lData.shape[1]) :
    lDataHist.SetBinContent(ic+1, ir+1, lData[ir,ic])

for ir in range(1, uDataHist.GetNbinsX()) :
  for ic in range(1, uDataHist.GetNbinsY()) :
    dNorm += abs(uDataHist.GetBinContent(ir, ic))

for ir in range(1, lDataHist.GetNbinsX()) :
  for ic in range(1, lDataHist.GetNbinsY()) :
    sNorm += abs(lDataHist.GetBinContent(ir, ic))
dNorm = 1300.0/dNorm
sNorm = 1300.0/sNorm

uDataHist.Scale(dNorm)
lDataHist.Scale(sNorm)

#Flip UED data sign, flipped back later
u1dData *= -1
u1dSclData *= -1
u1dSclSim *= -1

uDatMax = np.argmax(u1dSclSim)
uDatMin = np.argmin(u1dSclSim)
lDatMax = np.argmax(l1dSclSim)
lDatMin = np.argmin(l1dSclSim)
uDatRev = uDatMin + np.argmin(np.absolute(u1dSclSim[uDatMin:uDatMax] - u1dSclSim[-1]))
lDatRev = lDatMin + np.argmin(np.absolute(l1dSclSim[lDatMin:lDatMax] - l1dSclSim[-1]))


simAlnLng = 2000
uToff = 318.35
uTlow = 316.7
uThigh = 319.5
uDataHist.GetXaxis().SetLimits(uToff - (uDatRev+1)*0.1, uToff + (u1dData.shape[0]-(uDatRev+1)+1)*0.1)
uDataHist.SetContour(80)
u1dDataHist = ROOT.TH1F("UED1dData", "UED1dData", \
    u1dData.shape[0], uToff - (uDatRev+1)*0.1, uToff + (u1dData.shape[0]-(uDatRev+1)+1)*0.1)
u1dSclDataHist = ROOT.TH1F("UED1dSclData", "UED1dSclData", \
    u1dSclData.shape[0], uToff -(uDatRev+1)*0.1, uToff + (u1dSclData.shape[0]-(uDatRev+1)+1)*0.1)
u1dSclSimHist = ROOT.TH1F("UED1dSclSim", "UED1dSclSim", \
    u1dSclSim.shape[0], uToff - (uDatRev+1)*0.1, uToff + (u1dSclSim.shape[0]-(uDatRev+1)+1)*0.1)
u1dSimAlnHist = ROOT.TH1F("UED1dSimAln", "UED1dSimAln", simAlnLng, -simAlnLng/2000, simAlnLng/2000)
for it in range(u1dData.shape[0]) :
  u1dDataHist.SetBinContent(it+1, -1*u1dData[it])
u1dDataHist.SetError(u1dDataStd)
for it in range(u1dSclData.shape[0]) :
  u1dSclDataHist.SetBinContent(it+1, -1*u1dSclData[it])
for it in range(u1dSclSim.shape[0]) :
  u1dSclSimHist.SetBinContent(it+1, -1*u1dSclSim[it])
u1dSclSimHist.Scale(1.0/u1dDataHist.GetMaximum())
u1dDataHist.Scale(1.0/u1dDataHist.GetMaximum())

uSimAlnMax = np.argmax(u1dSimAln)
uSimAlnMin = np.argmin(u1dSimAln)
#uSimAlnRev = uSimAlnMin + np.argmin(np.absolute(u1dSimAln[uSimAlnMin:uSimAlnMax] - u1dSimAln[-1])) 
uSimAlnRev = np.argmin(np.absolute(u1dSimAln[0:uSimAlnMin] - 0.3333333)) 
for it in range(1, simAlnLng+1) :
  u1dSimAlnHist.SetBinContent(it, u1dSimAln[uSimAlnRev+it-simAlnLng/5.])


lToff = 318.72
l1dDataHist = ROOT.TH1F("LCLS1dData", "LCLS1dData", \
    l1dData.shape[0], lToff - (lDatRev+1)*4.936e-3, lToff + (l1dData.shape[0]-(lDatRev+1)+1)*4.936e-3)
l1dDataHist.SetError(l1dDataStd/np.sqrt(5.))
l1dSclDataHist = ROOT.TH1F("LCLS1dSclData", "LCLS1dSclData", \
    l1dSclData.shape[0], lToff - (lDatRev+1)*4.936e-3, lToff + (l1dSclData.shape[0]-(lDatRev+1)+1)*4.936e-3)
l1dSclSimHist = ROOT.TH1F("LCLS1dSclSim", "LCLS1dSclSim", \
    l1dSclSim.shape[0], lToff - (lDatRev+1)*4.936e-3, lToff + (l1dSclSim.shape[0]-(lDatRev+1)+1)*4.936e-3)
print("LCLS INFO: ",l1dSclData.shape[0],l1dSclData.shape[0],(lDatRev+1) + (l1dSclSim.shape[0]-(lDatRev+1)+1), lToff - (lDatRev+1)*4.936e-3,lToff + (l1dSclData.shape[0]-(lDatRev+1)+1)*4.936e-3)
sys.exit(0)
l1dSimAlnHist = ROOT.TH1F("LCLS1dSimAln", "LCLS1dSimAln", simAlnLng, -simAlnLng/2000, simAlnLng/2000)
l1dDataHist.GetXaxis().SetNdivisions(504)
l1dSclDataHist.GetXaxis().SetNdivisions(504)
l1dSclSimHist.GetXaxis().SetNdivisions(504)
for it in range(l1dData.shape[0]) :
  l1dDataHist.SetBinContent(it+1, l1dData[it])
for it in range(l1dSclData.shape[0]) :
  l1dSclDataHist.SetBinContent(it+1, l1dSclData[it])
for it in range(l1dSclSim.shape[0]) :
  l1dSclSimHist.SetBinContent(it+1, l1dSclSim[it])
l1dSclSimHist.Scale(1.0/l1dDataHist.GetMaximum())
l1dDataHist.Scale(1.0/l1dDataHist.GetMaximum())

lSimAlnMax = np.argmax(l1dSimAln)
lSimAlnMin = np.argmin(l1dSimAln)
lSimAlnRev = np.argmin(np.absolute(l1dSimAln[0:lSimAlnMin] - 0.33333333)) 
for it in range(1, simAlnLng+1) :
  l1dSimAlnHist.SetBinContent(it, l1dSimAln[lSimAlnRev+it-simAlnLng/5.])



####  Canvas ####
ypix = 500
xpix = 1200 
canReg = ROOT.TCanvas("canReg", "canReg", xpix, ypix);
hSep = 0.05
hrCut = (uThigh - uTlow)/((uThigh - uTlow) + (l1dSclSim.shape[0]-(lDatRev+1)+1)*4.936e-3 + (lDatRev+1)*4.936e-3)
ycut = 0.3
ltrScale = 0.85

pad1 = ROOT.TPad("pad1", "pad1", 0.0, ycut, hrCut*(1 - hSep), 1.0) 
pad1.SetBottomMargin(0.02)
pad1.SetLeftMargin(0.08)
pad1.SetRightMargin(0.1)
pad1.Draw()
pad1.cd()
labelSize = 0.05
#palette = uDataHist.GetListOfFunctions().FindObject("palette")
#palette.SetY1NDC(0.02)
#palette.SetX1NDC(0.835)
#palette.SetX2NDC(0.87)
norm = np.max(np.abs([uDataHist.GetMaximum(), uDataHist.GetMinimum()]))
uDataHist.GetYaxis().SetLabelSize(labelSize*1.15)
uDataHist.GetYaxis().SetTitleSize(labelSize*1.15)
uDataHist.GetYaxis().SetTitleOffset(0.6)
uDataHist.GetYaxis().SetTitle("Q[#AA^{-1}]")
uDataHist.GetZaxis().SetLabelSize(labelSize*1.15)
uDataHist.GetZaxis().SetTitleSize(labelSize*1.15)
uDataHist.Scale(norm)
#uDataHist.SetMaximum(uDataHist.GetMaximum()*norm)
#uDataHist.SetMinimum(uDataHist.GetMinimum()*norm)
plc.centerAxisTitles(uDataHist)
uDataHist.GetXaxis().SetRangeUser(uTlow, uThigh)
uDataHist.Draw("COLZ")

txt.SetTextSize(0.12*ltrScale)
txt.DrawText(319.4, 13.2, "a")


canReg.cd()
pad2 = ROOT.TPad("pad2", "pad2", hrCut*(1 + hSep), ycut, 1.0, 1.0) 
pad2.SetBottomMargin(0.02)
pad2.SetLeftMargin(0.0)
pad2.SetRightMargin(0.18)
pad2.Draw()
pad2.cd()
labelSize = 0.093

lDataHist.RebinX(5)
norm = 1./np.max(np.abs([lDataHist.GetMaximum(), lDataHist.GetMinimum()]))
lDataHist.GetYaxis().SetLabelSize(labelSize*0.85)
lDataHist.GetYaxis().SetTitleSize(labelSize*0.85)
lDataHist.GetYaxis().SetTitleOffset(1.5)
lDataHist.GetYaxis().SetTitle("Kinetic Energy [eV]")
lDataHist.GetZaxis().SetLabelSize(labelSize*0.6)
lDataHist.GetZaxis().SetTitleSize(labelSize*0.6)
lDataHist.Scale(norm)
plc.centerAxisTitles(lDataHist)
lDataHist.Draw("COLZ")

txt.SetTextSize(0.15*ltrScale)
#txt.DrawText(2.8, 528, "b")
txt.DrawText(2.85, 528, "b")


bMarg = 0.3
labelSize = 1.55*labelSize
canReg.cd()
pad3 = ROOT.TPad("pad3", "pad3", 0.0, 0, hrCut*(1 - hSep), ycut) 
pad3.Draw()
pad3.SetLeftMargin(0.08)
pad3.SetRightMargin(0.1)
pad3.SetBottomMargin(bMarg)
pad3.cd()
u1dDataHist.GetYaxis().SetNdivisions(505)
u1dDataHist.GetXaxis().SetLabelSize(labelSize)
u1dDataHist.GetXaxis().SetTitleSize(labelSize)
u1dDataHist.GetXaxis().SetTitle("Time [ps]")
u1dDataHist.GetYaxis().SetLabelSize(labelSize)
u1dDataHist.GetYaxis().SetTitleSize(labelSize)
u1dDataHist.GetYaxis().SetTitleOffset(0.22)
u1dDataHist.GetYaxis().SetTitle("P^{0}_{2} Coefficient")
u1dDataHist.SetMarkerStyle(20)
u1dDataHist.SetMarkerColor(4)
plc.centerAxisTitles(u1dDataHist)
legend3 = ROOT.TLegend(0.1,0.7,0.24,0.95)
legend3.AddEntry(u1dDataHist, "UED Data", "p")
legend3.AddEntry(u1dSclSimHist, "Fit: 55K", "l")
legend3.SetBorderSize(1)


u1dDataHist.GetXaxis().SetRangeUser(uTlow, uThigh)
u1dSclSimHist.GetXaxis().SetRangeUser(uTlow, uThigh)
u1dDataHist.Draw("PE")
u1dSclSimHist.Draw("SAMEC")
legend3.Draw()

txt.SetTextSize(0.27*ltrScale)
txt.DrawText(319.4, 0.84, "c")


canReg.cd()
pad4 = ROOT.TPad("pad4", "pad4", hrCut*(1 + hSep), 0, 1.0, ycut) 
pad4.Draw()
pad4.cd()
pad4.SetBottomMargin(bMarg)
pad4.SetLeftMargin(0.0)
pad4.SetRightMargin(0.18)
l1dDataHist.Rebin(5)
l1dDataHist.Scale(1.0/5.0)
l1dSclSimHist.Rebin(5)
l1dSclSimHist.Scale(1.0/5.0)
l1dDataHist.GetYaxis().SetNdivisions(505)
l1dDataHist.GetXaxis().SetLabelSize(labelSize)
l1dDataHist.GetXaxis().SetTitleSize(labelSize)
l1dDataHist.GetXaxis().SetTitle("Time [ps]")
l1dDataHist.GetYaxis().SetLabelSize(labelSize)
l1dDataHist.GetYaxis().SetTitleSize(labelSize)
l1dDataHist.GetYaxis().SetTitle("P^{0}_{2} Coefficient")
l1dDataHist.GetYaxis().SetTitleOffset(0.85)
l1dDataHist.SetMarkerStyle(20)
l1dDataHist.SetMarkerColor(2)
plc.centerAxisTitles(l1dDataHist)

legend4 = ROOT.TLegend(0.03,0.75,0.32,0.95)
legend4.AddEntry(l1dDataHist, "LCLS Data", "p")
legend4.AddEntry(l1dSclSimHist, "Fit: 300K", "l")
legend4.SetTextSize(0.075)
legend4.SetBorderSize(1)

l1dDataHist.Draw("PE")
l1dSclSimHist.Draw("SAMEC")
legend4.Draw()
pad4.RedrawAxis()

txt.SetTextSize(0.27*ltrScale)
txt.DrawText(lToff + 0.58, 0.8, "d")

"""
canReg.cd()
pad5 = ROOT.TPad("pad5", "pad5", 0.0, 0.03, 1.0, 0.205) 
pad5.Draw()
pad5.cd()
pad5.SetRightMargin(0.045)
pad5.SetLeftMargin(0.08)
labelSize = 0.14
u1dSimAlnHist.GetXaxis().SetLimits(-0.4, 1.6)
l1dSimAlnHist.GetXaxis().SetLimits(-0.4, 1.6)
u1dSimAlnHist.GetXaxis().SetLabelSize(labelSize)
u1dSimAlnHist.GetXaxis().SetTitleSize(labelSize)
u1dSimAlnHist.GetXaxis().SetTitle("Time [ps]")
u1dSimAlnHist.GetYaxis().SetLabelSize(labelSize)
u1dSimAlnHist.GetYaxis().SetTitleSize(labelSize)
u1dSimAlnHist.GetYaxis().SetTitleOffset(0.22)
u1dSimAlnHist.GetYaxis().SetTitle("\langle Cos^2(\\theta) \\rangle")
plc.centerAxisTitles(u1dSimAlnHist)

u1dSimAlnHist.SetLineColor(4)
l1dSimAlnHist.SetLineColor(2)

legend5 = ROOT.TLegend(0.1,0.61,0.24,0.91)
legend5.AddEntry(u1dSimAlnHist, "UED Fit: 55K", "l")
legend5.AddEntry(l1dSimAlnHist, "LCLS Fit: 300K", "l")
legend5.SetBorderSize(1)

u1dSimAlnHist.Draw()
l1dSimAlnHist.Draw("SAME")
legend5.Draw()

txt.SetTextSize(0.26)
txt.DrawText(1.5, 0.27, "E")
#txt.DrawText(1.5, 0.415, "E")


canReg.cd()
canReg.Update()
"""

canReg.Print("dataSimCompare.png")




#lSimAlnRev

plc = PLOTclass()
ROOT.gROOT.SetStyle("genStyle")
canComp = ROOT.TCanvas("canComp", "canComp", 800, 600);


rebin = 50
labelSize = 0.04
u1dSimAlnHist.Rebin(rebin)
l1dSimAlnHist.Rebin(rebin)
u1dSimAlnHist.Scale(1./rebin)
l1dSimAlnHist.Scale(1./rebin)

scale = 2.0
u1dSimAlnHist.GetXaxis().SetLimits(-simAlnLng/2000 - 0.4 + 0.025, simAlnLng/2000 - 0.4 - 0.025) 
l1dSimAlnHist.GetXaxis().SetLimits(-simAlnLng/2000*scale - 0.4 + 0.025, simAlnLng/2000*scale - 0.4 - 0.025) 
#u1dSimAlnHist.GetXaxis().SetLimits(-simAlnLng/2000 + (uSimAlnRev - uSimAlnMax)/1000, simAlnLng/2000 + (uSimAlnRev - uSimAlnMax)/1000) 
#l1dSimAlnHist

u1dSimAlnHist.SetLineWidth(5)
l1dSimAlnHist.SetLineWidth(5)
u1dSimAlnHist.GetXaxis().SetLabelSize(labelSize)
u1dSimAlnHist.GetXaxis().SetTitleSize(labelSize)
u1dSimAlnHist.GetXaxis().SetTitleOffset(1.4)
u1dSimAlnHist.GetXaxis().SetTitleColor(4)
u1dSimAlnHist.GetXaxis().SetTitle("Time [ps]")
u1dSimAlnHist.GetYaxis().SetLabelSize(labelSize)
u1dSimAlnHist.GetYaxis().SetTitleSize(labelSize)
u1dSimAlnHist.GetYaxis().SetTitleOffset(1.4)
u1dSimAlnHist.GetYaxis().SetTitle("\langle Cos^2(\\theta) \\rangle")
plc.centerAxisTitles(u1dSimAlnHist)

u1dSimAlnHist.GetXaxis().SetLabelColor(4)



legend5 = ROOT.TLegend(0.18,0.7,0.5,0.83)
legend5.AddEntry(u1dSimAlnHist, "UED Fit: 55K", "l")
legend5.AddEntry(l1dSimAlnHist, "LCLS Fit: 300K", "l")
legend5.SetTextSize(0.045)
legend5.SetBorderSize(1)

u1dSimAlnHist.Draw()
l1dSimAlnHist.Draw("SAME")
legend5.Draw()

"""
print(ROOT.gPad.GetUxmin(),ROOT.gPad.GetUymax(),ROOT.gPad.GetUxmax(),ROOT.gPad.GetUymax())
ax = ROOT.TGaxis(0.1, 0.1, 0.9, 0.1,1,2,50,"-")
ax.Draw()
"""

ROOT.gPad.SetTopMargin(0.13)

txt = ROOT.TText()
txt.SetTextSize(labelSize)
txt.SetTextAlign(22)
txt.SetTextColor(2)
axLabels = ["-0.6", "-0.5", "-0.4", "-0.3", "-0.2", "-0.1", "0 ", " 0.1", " 0.2"]
for i in range(len(axLabels)):
  #pos = (u1dSimAlnHist.GetBinCenter(i*4+4) + u1dSimAlnHist.GetBinCenter(i*5+4))/2.
  pos = u1dSimAlnHist.GetBinCenter(i*4+4) + 0.01 + 0.004*i 
  print(pos)
  txt.DrawText(pos, u1dSimAlnHist.GetMaximum()*1.037, axLabels[i])
txt.DrawText(-0.4, u1dSimAlnHist.GetMaximum()*1.077, "Time [ps]")



#ln = ROOT.TLine(0, u1dSimAlnHist.GetYaxis().GetXmin(), 0, u1dSimAlnHist.GetYaxis().GetXmax())
ln = ROOT.TLine(0, 0.232, 0, 0.452)
ln.SetLineWidth(4)
ln.Draw()

canComp.Print("SimBinCompare.png")

