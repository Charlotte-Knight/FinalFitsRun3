import ROOT
import numpy as np
import logging
log = logging.getLogger(__name__)

functions = ["Gaussian", "DCB", "Bernstein", "Exponential", "Power", "ExpPoly", "Chebychev", "Laurent"]

def setPrePostfix(*args, prefix="", postfix="", change_name=True, change_title=False):
  for arg in args:
    if change_name:
      arg.SetName(prefix+arg.GetName()+postfix)
    if change_title:
      arg.SetTitle(prefix+arg.GetTitle()+postfix)

class FinalFitsFunction:
  def __init__(self, x, prefix="", postfix="", bounds=None, order=1, randomise=False):
    self.order = order
    self.initBounds(bounds)
    self.initVars()
    self.initPdf(x)
    setPrePostfix(self.pdf, *self.vars, prefix=prefix, postfix=postfix)
    if randomise:
      self.randomiseVars()

  def initBounds(self, bounds):
    self.bounds = bounds if bounds is not None else self.getDefaultBounds() 
    
  def initVars(self):
    self.vars = [ROOT.RooRealVar(name, name, *self.bounds[name]) for name in self.bounds]
      
  def initPdf(self, x):
    self.pdf = self.getBasePdf(x, self.vars)
      
  def randomiseVars(self, skip_constant=True, seed=None):
    np.random.seed(seed)
    for v in self.vars:
      if (not v.isConstant()) or (not skip_constant):
        v.setVal(np.random.uniform(v.getMin(), v.getMax()))

  def getBasePdf(self, x, vars):
    pass

  def getDefaultBounds(self):
    pass
    
  def getDof(self):
    return len(self.vars)
  
  def checkBounds(self):
    for var in self.vars:
      if np.isclose(var.getVal(), var.getMin(), rtol=0.01) or np.isclose(var.getVal(), var.getMax(), rtol=0.01):
        log.warning(f"Parameter {var.GetName()} from function {self.pdf.GetName()} is at its bounds")
        log.warning(f"{var.GetName()}={var.getVal()}, low={var.getMin()}, high={var.getMax()}")

class FinalFitsFunctionSum(FinalFitsFunction):
  def initVars(self):
    if self.order == 1:
      super().initVars()
    else:
      self.vars = []

      self.params = [] # contains lists of parameters for each pdf
      for i in range(self.order):
        new_params = [ROOT.RooRealVar(name, name, *self.bounds[name]) for name in self.bounds]
        setPrePostfix(*new_params, postfix=f"{i}", change_title=True)
        self.params.append(new_params)
        self.vars += new_params

      self.coeffs = [ROOT.RooRealVar(f"c{i}", f"c{i}", 1/self.order, 0.0, 1.0) for i in range(self.order-1)]
      self.vars += self.coeffs
      
  def initPdf(self, x):
    if self.order == 1:
      super().initPdf(x)
    else:
      self.pdfs = [self.getBasePdf(x, params_subset) for params_subset in self.params]
      base_name = self.pdfs[0].GetName()
      base_title = self.pdfs[0].GetTitle()
      for i in range(self.order):
        setPrePostfix(self.pdfs[i], postfix=f"{i}")
      self.pdf = ROOT.RooAddPdf(f"{base_name}{self.order}", f"{base_title}{self.order}", self.pdfs, self.coeffs, True)
      self.pdf.fixCoefNormalization(x)

class Gaussian(FinalFitsFunctionSum):
  def getDefaultBounds(self):
    return {
      "mean": [125, 120, 130],
      "sigma": [1.5, 1, 5]
    }
  
  def getBasePdf(self, x, vars):
    return ROOT.RooGaussian("gauss", "Gaussian", x, *vars)

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
  
  def getBasePdf(self, x, vars):
    return ROOT.RooCrystalBall("dcb", "DCB", x, *vars)

class Bernstein(FinalFitsFunction):
  def getDefaultBounds(self):
    return {f"a{i}": [0.1, 0, 10] for i in range(self.order)}
  
  def getBasePdf(self, x, vars):
    return ROOT.RooBernstein(f"bern{self.order}", f"Bernstein{self.order}", x, [ROOT.RooFit.RooConst(1.0), *vars])
  
class Exponential(FinalFitsFunctionSum):
  def getDefaultBounds(self):
    return {"a": [-0.05, -0.1, 0]}
  
  def getBasePdf(self, x, vars):
    return ROOT.RooExponential("exp", "Exponential", x, *vars)

class Power(FinalFitsFunctionSum):
  def getDefaultBounds(self):
    return {"a": [-1, -5, 0]}

  def getBasePdf(self, x, vars):
    return ROOT.RooPower(f"pow", f"Power", x, *vars)

class ExpPoly(FinalFitsFunction):
  def getDefaultBounds(self):
    return {f"a{i}": [0.5, 0, 10] for i in range(self.order)}

  def getBasePdf(self, x, vars):
    formula = "+".join([f"@{i}*(@0/100)**{i}" for i in range(1, self.order+1)])
    self.poly = ROOT.RooFormulaVar("poly", "poly", formula, ROOT.RooArgList(x, *vars))
    return ROOT.RooExponential(f"exppoly{self.order}", f"ExpPoly{self.order}", self.poly, ROOT.RooFit.RooConst(-1.0))

class Chebychev(FinalFitsFunction):
  def getDefaultBounds(self):
    return {f"a{i}": [0.5, 0, 1] for i in range(self.order)}
  
  def getBasePdf(self, x, vars):
    return ROOT.RooChebychev(f"cheby{self.order}", f"Chebychev{self.order}", x, ROOT.RooArgList(*vars))
  
class Laurent(FinalFitsFunction):
  def getDefaultBounds(self):
    return {f"a{i}": [1/(self.order+1), 0, 1] for i in range(self.order)}
  
  def getBasePdf(self, x, vars):
    def g(i):
      return sum([(-1)**j * j for j in range(i+1)])
    
    self.pdfs = [ROOT.RooPower(f"pow{i}", f"Power{i}", x, ROOT.RooFit.RooConst(-4+g(i))) for i in range(self.order+1)]
    pdf = ROOT.RooAddPdf(f"laurent{self.order}", f"Laurent{self.order}", self.pdfs, vars, True)
    pdf.fixCoefNormalization(x)
    return pdf