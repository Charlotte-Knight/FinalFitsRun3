import ROOT

x = ROOT.RooRealVar("x", "x", -10, 10)
mean = ROOT.RooRealVar("mean", "mean", 0)
sigma = ROOT.RooRealVar("sigma", "sigma", 1)
gauss = ROOT.RooGaussian("gauss", "gauss", x, mean, sigma)