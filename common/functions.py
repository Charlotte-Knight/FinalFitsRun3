import ROOT
import numpy as np

def randomiseVars(vars, skip_constant=True, seed=1):
  np.random.seed(seed)
  for v in vars:
    if (not v.isConstant()) or (not skip_constant):
      v.setVal(np.random.uniform(v.getMin(), v.getMax()))

def setPrePostfix(*args, prefix="", postfix="", change_name=True, change_title=False):
  for arg in args:
    if change_name:
      arg.SetName(prefix+arg.GetName()+postfix)
    if change_title:
      arg.SetTitle(prefix+arg.GetTitle()+postfix)

class FinalFitsFunction:
  def __init__(self, x, prefix="", postfix="", bounds=None, order=1, randomise=False):
    self.initVars(bounds, randomise, order)
    self.initPdf(x, order)
    setPrePostfix(self.pdf, *self.vars, prefix=prefix, postfix=postfix)

  def initBounds(self, bounds):
    if bounds is None:
      self.bounds = self.getDefaultBounds()
    else:
      self.bounds = bounds
    
  def initVars(self, bounds, randomise, order):
    self.vars = []

    self.initBounds(bounds)
    self.params = [ROOT.RooRealVar(name, name, *self.bounds[name]) for name in self.bounds]
    self.vars += self.params
    if randomise:
      randomiseVars(self.vars)

  def getDefaultBounds(self):
    pass

  def getBasePdf(self, x, params, order):
    pass

  def initPdf(self, x, order):
    self.pdf = self.getBasePdf(x, self.params, order)

class FinalFitsFunctionSum(FinalFitsFunction):
  def initVars(self, bounds, randomise, order):
    if order == 1:
      super().initVars(bounds, randomise, order)
    else:
      self.vars = []
      self.initBounds(bounds)

      self.params = []
      for i in range(order):
        new_params = [ROOT.RooRealVar(name, name, *self.bounds[name]) for name in self.bounds]
        setPrePostfix(*new_params, postfix=f"{i}", change_title=True)
        self.params.append(new_params)
        self.vars += new_params

      self.coeffs = [ROOT.RooRealVar(f"c{i}", f"c{i}", 1/order, 0.00, 1.00) for i in range(order-1)]
      self.vars += self.coeffs

      if randomise:
        randomiseVars(self.vars)
      
  def initPdf(self, x, order):
    if order == 1:
      super().initPdf(x, order)
    else:
      self.pdfs = [self.getBasePdf(x, params_subset, order=1) for params_subset in self.params]
      for i in range(order):
        setPrePostfix(self.pdfs[i], postfix=f"{i}")
      base_name = self.pdfs[0].GetName()[:-1]
      base_title = self.pdfs[0].GetTitle()
      self.pdf = ROOT.RooAddPdf(f"{base_name}{order}", f"{base_title}{order}", ROOT.RooArgList(*self.pdfs), ROOT.RooArgList(*self.coeffs), True)
      self.pdf.fixCoefNormalization(ROOT.RooArgSet(x))

class Gaussian(FinalFitsFunctionSum):
  def getDefaultBounds(self):
    return {
      "mean": [125, 120, 130],
      "sigma": [1.5, 1, 5]
    }
  
  def getBasePdf(self, x, params, order):
    return ROOT.RooGaussian("gauss", "Gaussian", x, *params)

class DCB(FinalFitsFunction):
  def getDefaultBounds(self):
    return {
      "mean": [125, 120, 130],
      "sigmaLR": [1.3, 1, 5],
      "alphaL": [1, 0, 5],
      "nL": [5, 0, 20],
      "alphaR": [1, 0, 5],
      "nR": [20, 0, 20]
    }
  
  def getBasePdf(self, x, params, order):
    return ROOT.RooCrystalBall("dcb", "DCB", x, *params)

class Bernstein(FinalFitsFunction):
  def getDefaultBounds(self):
    return {
      "a": [0.1, 0, 10]
    }
  
  def initVars(self, bounds, randomise, order):
    self.vars = []

    self.initBounds(bounds)
    self.params = [ROOT.RooRealVar(f"a{i}", f"a{i}", *self.bounds["a"]) for i in range(order)]
    self.vars += self.params
    if randomise:
      randomiseVars(self.vars)

  def getBasePdf(self, x, params, order):
    return ROOT.RooBernstein(f"bern{order}", f"Bernstein{order}", x, ROOT.RooArgList(ROOT.RooFit.RooConst(1.0), *params))

# class BernsteinFast(FinalFitsFunction):
#   def getDefaultBounds(self):
#     return {
#       "a": [0.1, 0, 10]
#     }
  
#   def initVars(self, bounds, randomise, order):
#     self.vars = []

#     self.initBounds(bounds)
#     self.params = [ROOT.RooRealVar(f"a{i}", f"a{i}", *self.bounds["a"]) for i in range(order)]
#     self.vars += self.params
#     if randomise:
#       randomiseVars(self.vars)

#   def getBasePdf(self, x, params, order):
#     return ROOT.RooBernsteinFast(f"bern{order}", f"Bernstein{order}", x, ROOT.RooArgList(*params))


class Exponential(FinalFitsFunctionSum):
  def getDefaultBounds(self):
    return {
      "a": [-0.05, -0.1, 0],
    }
  
  def getBasePdf(self, x, params, order):
    return ROOT.RooExponential("exp", "Exponential", x, *params)

class Power(FinalFitsFunctionSum):
  def getDefaultBounds(self):
    return {
      "a": [-1, -10, 0],
    }

  def getBasePdf(self, x, params, order):
    return ROOT.RooPower(f"pow", f"Power", x, *params)

# class ExpPoly(FinalFitsFunction):
#   def getDefaultBounds(self):
#     return {
#       "a": [-0.1, -1, -0.00001],
#     }

#   def initVars(self, bounds, randomise, order):
#     self.vars = []

#     self.initBounds(bounds)
#     self.params = [ROOT.RooRealVar(f"a{i}", f"a{i}", *self.bounds["a"]) for i in range(order)]
#     self.vars += self.params
#     if randomise:
#       randomiseVars(self.vars)

#   def getBasePdf(self, x, params, order):
#     return ROOT.RooExpPoly(f"exppoly{order}", f"ExpPoly{order}", x, ROOT.RooArgList(*params))

class ExpPoly(FinalFitsFunction):
  def getDefaultBounds(self):
    return {
      "a": [0.5, 0, 2],
    }

  def initVars(self, bounds, randomise, order):
    self.vars = []

    self.initBounds(bounds)
    self.params = [ROOT.RooRealVar(f"a{i}", f"a{i}", *self.bounds["a"]) for i in range(order)]
    self.vars += self.params
    if randomise:
      randomiseVars(self.vars)

  def getBasePdf(self, x, params, order):
    formula = "+".join([f"@{i+1}*(@0/100)" for i in range(order)])
    self.poly = ROOT.RooFormulaVar("poly", "poly", formula, ROOT.RooArgList(x, *params))
    return ROOT.RooExponential(f"exppoly{order}", f"ExpPoly{order}", self.poly, ROOT.RooFit.RooConst(-1.0))

# class Laurent(FinalFitsFunction):
#   def getDefaultBounds(self):
#     return {

#     }

class Chebychev(FinalFitsFunction):
  def getDefaultBounds(self):
    return {
      "a": [0.5, 0, 1]
    }
  
  def initVars(self, bounds, randomise, order):
    self.vars = []

    self.initBounds(bounds)
    self.params = [ROOT.RooRealVar(f"a{i}", f"a{i}", *self.bounds["a"]) for i in range(order)]
    self.vars += self.params
    if randomise:
      randomiseVars(self.vars)

  def getBasePdf(self, x, params, order):
    return ROOT.RooChebychev(f"cheby{order}", f"Chebychev{order}", x, ROOT.RooArgList(*params))