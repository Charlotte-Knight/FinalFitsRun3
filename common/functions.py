import ROOT
import numpy as np
import logging
log = logging.getLogger(__name__)

functions = ["Gaussian", "DCB", "Bernstein", "Exponential", "Power", "ExpPoly", "Chebychev"]

def randomiseVars(vars, skip_constant=True, seed=None):
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
    self.order = order
    self.initVars(bounds, randomise)
    self.initPdf(x)
    setPrePostfix(self.pdf, *self.vars, prefix=prefix, postfix=postfix)

  def initBounds(self, bounds):
    if bounds is None:
      self.bounds = self.getDefaultBounds()
    else:
      self.bounds = bounds
    
  def initVars(self, bounds, randomise):
    self.vars = []

    self.initBounds(bounds)
    self.params = [ROOT.RooRealVar(name, name, *self.bounds[name]) for name in self.bounds]
    self.vars += self.params
    if randomise:
      randomiseVars(self.vars)

  def getDefaultBounds(self):
    pass

  def getBasePdf(self, x, params):
    pass

  def initPdf(self, x):
    self.pdf = self.getBasePdf(x, self.params)
  
  def getDof(self):
    return len(self.vars)

class FinalFitsFunctionSum(FinalFitsFunction):
  def initVars(self, bounds, randomise):
    if self.order == 1:
      super().initVars(bounds, randomise)
    else:
      self.vars = []
      self.initBounds(bounds)

      self.params = []
      for i in range(self.order):
        new_params = [ROOT.RooRealVar(name, name, *self.bounds[name]) for name in self.bounds]
        setPrePostfix(*new_params, postfix=f"{i}", change_title=True)
        self.params.append(new_params)
        self.vars += new_params

      self.coeffs = [ROOT.RooRealVar(f"c{i}", f"c{i}", 1/self.order, 0.0, 1.0) for i in range(self.order-1)]
      self.vars += self.coeffs

      if randomise:
        randomiseVars(self.vars)
      
  def initPdf(self, x):
    if self.order == 1:
      super().initPdf(x)
    else:
      self.pdfs = [self.getBasePdf(x, params_subset) for params_subset in self.params]
      for i in range(self.order):
        setPrePostfix(self.pdfs[i], postfix=f"{i}")
      base_name = self.pdfs[0].GetName()[:-1]
      base_title = self.pdfs[0].GetTitle()
      self.pdf = ROOT.RooAddPdf(f"{base_name}{self.order}", f"{base_title}{self.order}", ROOT.RooArgList(*self.pdfs), ROOT.RooArgList(*self.coeffs), True)
      self.pdf.fixCoefNormalization(ROOT.RooArgSet(x))

class Gaussian(FinalFitsFunctionSum):
  def getDefaultBounds(self):
    return {
      "mean": [125, 120, 130],
      "sigma": [1.5, 1, 5]
    }
  
  def getBasePdf(self, x, params):
    return ROOT.RooGaussian("gauss", "Gaussian", x, *params)

class DCB(FinalFitsFunction):
  def getDefaultBounds(self):
    return {
      "mean": [125, 120, 130],
      "sigmaLR": [1.3, 1, 5],
      "alphaL": [1, 0.1, 5],
      "nL": [5, 0.1, 20],
      "alphaR": [1, 0.1, 5],
      "nR": [20, 0.1, 20]
    }
  
  def getBasePdf(self, x, params):
    return ROOT.RooCrystalBall("dcb", "DCB", x, *params)

class Bernstein(FinalFitsFunction):
  def getDefaultBounds(self):
    return {
      "a": [0.1, 0, 10]
    }
  
  def initVars(self, bounds, randomise):
    self.vars = []

    self.initBounds(bounds)
    self.params = [ROOT.RooRealVar(f"a{i}", f"a{i}", *self.bounds["a"]) for i in range(self.order)]
    self.vars += self.params
    if randomise:
      randomiseVars(self.vars)

  def getBasePdf(self, x, params):
    return ROOT.RooBernstein(f"bern{self.order}", f"Bernstein{self.order}", x, ROOT.RooArgList(ROOT.RooFit.RooConst(1.0), *params))
  
class BernsteinFast(FinalFitsFunction):
  def getDefaultBounds(self):
    return {
      "a": [0.1, 0, 10]
    }
  
  def initVars(self, bounds, randomise):
    self.vars = []

    self.initBounds(bounds)
    self.params = [ROOT.RooRealVar(f"a{i}", f"a{i}", *self.bounds["a"]) for i in range(self.order)]
    self.vars += self.params
    if randomise:
      randomiseVars(self.vars)

  def getBasePdf(self, x, params):
    return ROOT.RooBernsteinFast(self.order)(f"bern{self.order}", f"Bernstein{self.order}", x, ROOT.RooArgList(*params))

class Exponential(FinalFitsFunctionSum):
  def getDefaultBounds(self):
    return {
      "a": [-0.05, -0.1, 0],
    }
  
  def getBasePdf(self, x, params):
    return ROOT.RooExponential("exp", "Exponential", x, *params)

class Power(FinalFitsFunctionSum):
  def getDefaultBounds(self):
    return {
      "a": [-1, -5, 0],
    }

  def getBasePdf(self, x, params):
    return ROOT.RooPower(f"pow", f"Power", x, *params)

class ExpPoly(FinalFitsFunction):
  def getDefaultBounds(self):
    return {
      "a": [0.5, 0, 10],
    }

  def initVars(self, bounds, randomise):
    self.vars = []

    self.initBounds(bounds)
    self.params = [ROOT.RooRealVar(f"a{i}", f"a{i}", *self.bounds["a"]) for i in range(self.order)]
    self.vars += self.params
    if randomise:
      randomiseVars(self.vars)

  def getBasePdf(self, x, params):
    formula = "+".join([f"@{i}*(@0/100)**{i}" for i in range(1, self.order+1)])
    self.poly = ROOT.RooFormulaVar("poly", "poly", formula, ROOT.RooArgList(x, *params))
    return ROOT.RooExponential(f"exppoly{self.order}", f"ExpPoly{self.order}", self.poly, ROOT.RooFit.RooConst(-1.0))

# class Laurent(FinalFitsFunction):
#   def getDefaultBounds(self):
#     return {

#     }

class Chebychev(FinalFitsFunction):
  def getDefaultBounds(self):
    return {
      "a": [0.5, 0, 1]
    }
  
  def initVars(self, bounds, randomise):
    self.vars = []

    self.initBounds(bounds)
    self.params = [ROOT.RooRealVar(f"a{i}", f"a{i}", *self.bounds["a"]) for i in range(self.order)]
    self.vars += self.params
    if randomise:
      randomiseVars(self.vars)

  def getBasePdf(self, x, params):
    return ROOT.RooChebychev(f"cheby{self.order}", f"Chebychev{self.order}", x, ROOT.RooArgList(*params))