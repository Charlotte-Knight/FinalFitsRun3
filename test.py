# import ROOT

# x = ROOT.RooRealVar("x", "x", 0, -10, 10)
# m = ROOT.RooRealVar("m", "m", 1, -10, 10)
# s = ROOT.RooRealVar("s", "s", 1, 0, 10)

# g = ROOT.RooGaussian("g", "g", x, m, s)

# x_c = x.Clone("x_c")
# data = g.generate({x_c}, 1000)
# data.SetName("data")

# w = ROOT.RooWorkspace("w", "w")
# w.Import(g)

# new_w = ROOT.RooWorkspace("new_w", "new_w")
# new_w.Import(w.pdf("g"))
# new_w.Import(data, ROOT.RooFit.RecycleConflictNodes())

# res = new_w.pdf("g").fitTo(new_w.data("data"), ROOT.RooFit.Save(), PrintLevel=-1)
# res.Print()

# w.pdf("g").Print("t")
# new_w.pdf("g").Print("t")

# w.Import(res)

import ROOT
from finalfits import pdfs

x = ROOT.RooRealVar("x", "x", 110, 140)
MH = ROOT.RooRealVar("MH", "MH", 125)
MH.setConstant(True)
pdf = getattr(pdfs, "Gaussian")(x, order=2,
                              transforms={"mean*": [MH, 1]}, polys={"mean*": [MH, 1]})