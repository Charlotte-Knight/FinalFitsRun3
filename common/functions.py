import ROOT
import numpy as np

def handlePrePostfix(*args, prefix="", postfix="", change_name=True, change_title=False):
  for arg in args:
    if change_name:
      arg.SetName(prefix+arg.GetName()+postfix)
    if change_title:
      arg.SetTitle(prefix+arg.GetTitle()+postfix)

def Gaussian(x, prefix="", postfix="", bounds=None, order=1):
  if bounds is None:
    bounds = {
      "mean": [125, 120, 130],
      "sigma": [1.5, 1, 5]
    }
  params = [ROOT.RooRealVar(name, name, *bounds[name]) for name in bounds]
  pdf = ROOT.RooGaussian("gauss", "Gaussian", x, *params)
  handlePrePostfix(pdf, *params, prefix=prefix, postfix=postfix)
  return pdf, params

def SumGaussians(x, prefix="", postfix="", bounds=None, order=2):
  if bounds is None:
    bounds = {
      "mean": [125, 120, 130],
      "sigma": [1.5, 0.5, 10]
    }
  all_params = []
  pdfs = []
  for i in range(order):
    pdf, params = Gaussian(x, bounds=bounds)
    handlePrePostfix(pdf, *params, postfix=f"{i}", change_title=True)
    all_params += params
    pdfs.append(pdf)

  for p in all_params:
    p.setVal(np.random.uniform(p.getMin(), p.getMax()))

  coeffs = [ROOT.RooRealVar(f"c{i}", f"c{i}", 1/order, -1.0,1.0) for i in range(order-1)]
  
  pdf = ROOT.RooAddPdf(f"gauss{order}", f"Gaussian{order}", ROOT.RooArgList(*pdfs), ROOT.RooArgList(*coeffs), True)
  pdf.fixCoefNormalization(ROOT.RooArgSet(x))

  other_objects = pdfs + all_params + coeffs
  handlePrePostfix(pdf, *other_objects, prefix=prefix, postfix=postfix)
  return pdf, other_objects


def DCB(x, prefix="", postfix="", bounds=None, order=1):
  if bounds is None:
    bounds = {
      "x0": [125, 120, 130],
      "sigmaLR": [1.3, 1, 5],
      "alphaL": [1, 0, 5],
      "nL": [5, 0, 20],
      "alphaR": [1, 0, 5],
      "nR": [20, 0, 20]
    }
  params = [ROOT.RooRealVar(name, name, *bounds[name]) for name in bounds]
  pdf = ROOT.RooCrystalBall("dcb", "DCB", x, *params)
  handlePrePostfix(pdf, *params, prefix=prefix, postfix=postfix)
  return pdf, params
  
  