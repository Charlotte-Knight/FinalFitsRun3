import ROOT
ROOT.gROOT.SetBatch(True)

 
# Create model for physics sample
# -------------------------------------------------------------
 
# Create observables
print("defining vars")
x = ROOT.RooRealVar("x", "x", 110, 140)
 
# Construct signal pdf
sigma = ROOT.RooRealVar("sigma", "sigma", 1, 0.5, 1.5)

a = ROOT.RooRealVar("a", "a", -0.1, -20, 20)
b = ROOT.RooRealVar("b", "b", 0.01, -1, 1)

MH1 = ROOT.RooFit.RooConst(120)
MH2 = ROOT.RooFit.RooConst(125)
MH3 = ROOT.RooFit.RooConst(130)

dm1 = ROOT.RooFormulaVar("dm1", "dm1", "@0 + @1*(@2-@3)", [a,b,MH1,MH2])
dm2 = ROOT.RooFormulaVar("dm2", "dm2", "@0 + @1*(@2-@3)", [a,b,MH2,MH2])
dm3 = ROOT.RooFormulaVar("dm3", "dm3", "@0 + @1*(@2-@3)", [a,b,MH3,MH2])

mean1 = ROOT.RooFormulaVar("mean1", "mean1", "@0 + @1", [MH1, dm1])
mean2 = ROOT.RooFormulaVar("mean2", "mean2", "@0 + @1", [MH2, dm2])
mean3 = ROOT.RooFormulaVar("mean3", "mean3", "@0 + @1", [MH3, dm3])

print(mean1)
print(mean2)
print(mean3)

print("init pdfs")
import copy
g1 = ROOT.RooGaussian("g1", "g1", x, mean1, sigma)
g2 = ROOT.RooGaussian("g2", "g2", x, mean2, sigma)
g3 = ROOT.RooGaussian("g3", "g3", x, mean3, sigma)
 
 
# Generate events for both samples
# ---------------------------------------------------------------
 
# Generate 1000 events in x and y from model
print("generating data")
# data1 = g1.generate({x}, 100)
# data2 = g2.generate({x}, 100)
# data3 = g3.generate({x}, 100)
data1 = g1.generateBinned({x}, 100)
data2 = g2.generateBinned({x}, 100)
data3 = g3.generateBinned({x}, 100)
 
# Create index category and join samples
# ---------------------------------------------------------------------------
 
# Define category to distinguish physics and control samples events
print("defining categories")
sample = ROOT.RooCategory("sample", "sample")
sample.defineType("one")
sample.defineType("two")
sample.defineType("three")
 
# Construct combined dataset in (x,sample)
print("combining data")
# combData = ROOT.RooDataSet(
#     "combData",
#     "combined data",
#     {x},
#     Index=sample,
#     Import={"one": data1, "two": data2, "three": data3}
# )
combData = ROOT.RooDataHist(
    "combData",
    "combined data",
    ROOT.RooArgList(x),
    Index=sample,
    Import={"one": data1, "two": data2, "three": data3}
)
 
# Construct a simultaneous pdf in (x, sample)
# -----------------------------------------------------------------------------------
 
# Construct a simultaneous pdf using category sample as index: associate model
# with the physics state and model_ctl with the control state
print("init sim pdf")
simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", {"one": g1, "two": g2, "three": g3}, sample)
 
# Perform a simultaneous fit
# ---------------------------------------------------

b.setConstant(True)
b.setVal(0)

dms = []
dms_err = []

r1 = g1.fitTo(data1, PrintLevel=-1, Save=True)
dms.append(dm1.getVal())
dms_err.append(a.getError())

r2 = g2.fitTo(data2, PrintLevel=-1, Save=True)
dms.append(dm2.getVal())
dms_err.append(a.getError())

r3 = g3.fitTo(data3, PrintLevel=-1, Save=True)
dms.append(dm3.getVal())
dms_err.append(a.getError())



b.setConstant(False)

# Perform simultaneous fit of model to data and model_ctl to data_ctl
print("do sim fit")
fitResult = simPdf.fitTo(combData, PrintLevel=-1, Save=True)
fitResult.Print()
 
# Plot model slices on data slices
# ----------------------------------------------------------------
 
# Make a frame for the physics sample
frame1 = x.frame(Title="120")
combData.plotOn(frame1, Cut="sample==sample::one")
simPdf.plotOn(frame1, Slice=(sample, "one"), ProjWData=(sample, combData))

frame2 = x.frame(Title="125")
combData.plotOn(frame2, Cut="sample==sample::two")
simPdf.plotOn(frame2, Slice=(sample, "two"), ProjWData=(sample, combData))

frame3 = x.frame(Title="130")
combData.plotOn(frame3, Cut="sample==sample::three")
simPdf.plotOn(frame3, Slice=(sample, "three"), ProjWData=(sample, combData))

c = ROOT.TCanvas("rf501_simultaneouspdf", "rf501_simultaneouspdf", 1200, 400)
c.Divide(3)
 
 
def draw(i, frame):
    c.cd(i)
    ROOT.gPad.SetLeftMargin(0.15)
    frame.GetYaxis().SetTitleOffset(1.4)
    frame.Draw()
 
 
draw(1, frame1)
draw(2, frame2)
draw(3, frame3)
 
c.SaveAs("rf501_simultaneouspdf.png")
 
import numpy as np
import matplotlib.pyplot as plt

print(dms)
print(dms_err)
xi = np.linspace(115, 135, 100)
plt.errorbar([MH1.getVal(), MH2.getVal(), MH3.getVal()], dms, dms_err, fmt=".", label="data")
plt.plot(xi, a.getVal() + b.getVal()*(xi-MH2.getVal()), label="simul fit")

p = np.polyfit([MH1.getVal(), MH2.getVal(), MH3.getVal()], dms, 1)
plt.plot(xi, p[1]+p[0]*xi, label="poly fit")
plt.legend()

plt.savefig("fit.png")

r1.Print()
r2.Print()
r3.Print()