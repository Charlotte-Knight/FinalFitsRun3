"""Classes that wrap around RooFit pdfs and their parameters."""

from typing import Optional
import logging
import ROOT
import numpy as np

log = logging.getLogger(__name__)

functions = ["Gaussian", "DCB", "Bernstein", "Exponential", "Power", "ExpPoly", "Chebychev",
             "Laurent"]

def set_pre_postfix(*args: tuple[ROOT.TNamed, ...], prefix: str = "", postfix: str = "",
                    change_name: bool = True, change_title: bool = False) -> None:
  """Set prefix and postfix for names and titles of RooFit objects
  
  Args:
      *args (Instances of ROOT TNamed objects): objects to change name/title of
      prefix (str, optional): string to append to start of object(s) name/title. Defaults to "".
      postfix (str, optional): string to append to end of object(s) name/title. Defaults to "".
      change_name (bool, optional): change object(s) name. Defaults to True.
      change_title (bool, optional): change object(s) title. Defaults to False.
  """
  for arg in args:
    if change_name:
      arg.SetName(prefix+arg.GetName()+postfix)
    if change_title:
      arg.SetTitle(prefix+arg.GetTitle()+postfix)

class FinalFitsPdf:
  """Base class for pdfs. Wraps around RooAbsPdf objects and their parameters."""
  def __init__(self, x: ROOT.RooRealVar, prefix: str = "", postfix: str = "",
               bounds: Optional[dict[str, tuple[float, float, float]]] = None, order: int = 1,
               randomise: bool = False) -> None:
    """
    Args:
        x (ROOT.RooRealVar): the random variable (usually mass) of the function
        prefix (str, optional): string to append to start of all RooFit object(s) name/title. Defaults to "".
        postfix (str, optional): string to append to start of all RooFit object(s) name/title. Defaults to "".
        bounds (Optional[dict[str, tuple[float, float, float]]], optional): Parameter bounds. Dictionary with parameter names as keys and tuples with (default value, min, max) as values. Defaults to None.
        order (int, optional): Order of the pdf. Defaults to 1.
        randomise (bool, optional): Randomise the values of the pdf parameters. Defaults to False.
    """
    self.order = order
    self.init_parameter_bounds(bounds)
    self.init_params()
    self.init_pdf(x)
    set_pre_postfix(self.pdf, *self.params, prefix=prefix, postfix=postfix)
    if randomise:
      self.randomise_params()

  def init_parameter_bounds(self, bounds: dict[str, tuple[float, float, float]]) -> None:
    """Initialize bounds for parameters of the pdf

    Args:
        bounds (dict[str, tuple[float, float, float]]): Dictionary with parameter names as keys and tuples with (default value, min, max) as values.
    """
    self.bounds = bounds if bounds is not None else self.get_default_bounds()

  def init_params(self):
    """Initialize parameters (RooRealVars) of the pdf"""
    self.params = [ROOT.RooRealVar(name, name, *self.bounds[name]) for name in self.bounds]

  def init_pdf(self, x: ROOT.RooRealVar) -> None:
    """Initialize the RooAbsPdf object

    Args:
        x (ROOT.RooRealVar): the random variable (usually mass) of the function
    """
    self.pdf = self.get_base_pdf(x, self.params)

  def randomise_params(self, skip_constant: bool = True, seed: bool = None) -> None:
    """Randomise the parameters of the pdf

    Args:
        skip_constant (bool, optional): skip variables which are set to be constant in RooFit. Defaults to True.
        seed (bool, optional): random seed. Defaults to None.
    """
    np.random.seed(seed)
    for v in self.params:
      if (not v.isConstant()) or (not skip_constant):
        v.setVal(np.random.uniform(v.getMin(), v.getMax()))

  def get_base_pdf(self, x: ROOT.RooRealVar, params: list[ROOT.RooRealVar]) -> ROOT.RooAbsPdf:
    """Get the "base" pdf object. For simple pdfs, this is the pdf itself. For more complex pdfs, such as a sum of pdfs, this is a component pdf. 

    Args:
        x (ROOT.RooRealVar): the random variable (usually mass) of the pdf
        params (list[ROOT.RooRealVar]): the parameters of the pdf

    Returns:
        ROOT.RooAbsPdf: base pdf object
    """
    pass

  def get_default_bounds(self) -> dict[str, tuple[float, float, float]]:
    """Get the default bounds for the parameters of the pdf

    Returns:
        dict[str, tuple[float, float, float]]: Dictionary with parameter names as keys and tuples with (default value, min, max) as values.
    """
    pass

  def get_dof(self) -> int:
    """Get the degrees of freedom of the pdf. This is the number of parameters of the pdf.

    Returns:
        int: degrees of freedom of the pdf
    """
    return len(self.params)

  def check_bounds(self) -> None:
    """Check if any of the parameters are at their bounds and log a warning if they are."""
    for var in self.params:
      if np.isclose(var.getVal(), var.getMin(), rtol=0.01) or np.isclose(var.getVal(), var.getMax(), rtol=0.01):
        log.warning("Parameter %s from function %s is at its bounds", var.GetName(), self.pdf.GetName())
        log.warning("%s=%s, low=%f, high=%f", var.GetName(), var.getVal(), var.getMin(), var.getMax())

class FinalFitsPdfSum(FinalFitsPdf):
  """Base class for pdfs that are sums of other pdfs."""
  def init_params(self) -> None:
    if self.order == 1:
      super().init_params()
    else:
      self.params = []

      # contains params for each pdf in sum w/o the coeffs that encode the relative contribitions
      self.params_no_coef = []
      for i in range(self.order):
        new_params = [ROOT.RooRealVar(name, name, *self.bounds[name]) for name in self.bounds]
        set_pre_postfix(*new_params, postfix=f"{i}", change_title=True)
        self.params_no_coef.append(new_params)
        self.params += new_params

      self.coeffs = [ROOT.RooRealVar(f"c{i}", f"c{i}", 1/self.order, 0.0, 1.0) for i in range(self.order-1)]
      self.params += self.coeffs

  def init_pdf(self, x: ROOT.RooRealVar) -> None:
    if self.order == 1:
      super().init_pdf(x)
    else:
      self.pdfs = [self.get_base_pdf(x, params_subset) for params_subset in self.params_no_coef]
      base_name = self.pdfs[0].GetName()
      base_title = self.pdfs[0].GetTitle()
      for i in range(self.order):
        set_pre_postfix(self.pdfs[i], postfix=f"{i}")
      self.pdf = ROOT.RooAddPdf(f"{base_name}{self.order}", f"{base_title}{self.order}", self.pdfs, self.coeffs, True)
      self.pdf.fixCoefNormalization(x)

class Gaussian(FinalFitsPdfSum):
  """Sum of Gaussians pdf"""
  def get_default_bounds(self) -> dict[str, tuple[float, float, float]]:
    return {
      "mean": [125, 120, 130],
      "sigma": [1.5, 1, 5]
    }

  def get_base_pdf(self, x: ROOT.RooRealVar, params: list[ROOT.RooRealVar]) -> ROOT.RooAbsPdf:
    return ROOT.RooGaussian("gauss", "Gaussian", x, *params)

class DCB(FinalFitsPdf):
  """Double Crystal Ball pdf"""
  def get_default_bounds(self) -> dict[str, tuple[float, float, float]]:
    return {
      "mean": [125, 120, 130],
      "sigmaLR": [1.3, 1, 5],
      "alphaL": [1, 0.1, 5],
      "nL": [5, 0.1, 20],
      "alphaR": [1, 0.1, 5],
      "nR": [20, 0.1, 20]
    }

  def get_base_pdf(self, x: ROOT.RooRealVar, params: list[ROOT.RooRealVar]) -> ROOT.RooAbsPdf:
    return ROOT.RooCrystalBall("dcb", "DCB", x, *params)

class Bernstein(FinalFitsPdf):
  """Bernstein polynomial pdf"""
  def get_default_bounds(self) -> dict[str, tuple[float, float, float]]:
    return {f"a{i}": [0.1, 0, 10] for i in range(self.order)}

  def get_base_pdf(self, x: ROOT.RooRealVar, params: list[ROOT.RooRealVar]) -> ROOT.RooAbsPdf:
    return ROOT.RooBernstein(f"bern{self.order}", f"Bernstein{self.order}", x, [ROOT.RooFit.RooConst(1.0), *params])

class Exponential(FinalFitsPdfSum):
  """Sum of exponentials pdf"""
  def get_default_bounds(self) -> dict[str, tuple[float, float, float]]:
    return {"a": [-0.05, -0.1, 0]}

  def get_base_pdf(self, x: ROOT.RooRealVar, params: list[ROOT.RooRealVar]) -> ROOT.RooAbsPdf:
    return ROOT.RooExponential("exp", "Exponential", x, *params)

class Power(FinalFitsPdfSum):
  """Sum of power laws pdf"""
  def get_default_bounds(self) -> dict[str, tuple[float, float, float]]:
    return {"a": [-1, -5, 0]}

  def get_base_pdf(self, x: ROOT.RooRealVar, params: list[ROOT.RooRealVar]) -> ROOT.RooAbsPdf:
    return ROOT.RooPower("pow", "Power", x, *params)

class ExpPoly(FinalFitsPdf):
  """Exponential polynomial pdf"""
  def get_default_bounds(self) -> dict[str, tuple[float, float, float]]:
    return {f"a{i}": [0.5, 0, 10] for i in range(self.order)}

  def get_base_pdf(self, x: ROOT.RooRealVar, params: list[ROOT.RooRealVar]) -> ROOT.RooAbsPdf:
    formula = "+".join([f"@{i}*(@0/100)**{i}" for i in range(1, self.order+1)])
    self.poly = ROOT.RooFormulaVar("poly", "poly", formula, ROOT.RooArgList(x, *params))
    return ROOT.RooExponential(f"exppoly{self.order}", f"ExpPoly{self.order}", self.poly, ROOT.RooFit.RooConst(-1.0))

class Chebychev(FinalFitsPdf):
  """Chebychev polynomial pdf"""
  def get_default_bounds(self) -> dict[str, tuple[float, float, float]]:
    return {f"a{i}": [0.5, 0, 1] for i in range(self.order)}

  def get_base_pdf(self, x: ROOT.RooRealVar, params: list[ROOT.RooRealVar]) -> ROOT.RooAbsPdf:
    return ROOT.RooChebychev(f"cheby{self.order}", f"Chebychev{self.order}", x, ROOT.RooArgList(*params))

class Laurent(FinalFitsPdf):
  """Laurent polynomial pdf"""
  def get_default_bounds(self) -> dict[str, tuple[float, float, float]]:
    return {f"a{i}": [1/(self.order+1), 0, 1] for i in range(self.order)}

  def get_base_pdf(self, x: ROOT.RooRealVar, params: list[ROOT.RooRealVar]) -> ROOT.RooAbsPdf:
    def g(i):
      return sum([(-1)**j * j for j in range(i+1)])

    self.pdfs = [ROOT.RooPower(f"pow{i}", f"Power{i}", x, ROOT.RooFit.RooConst(-4+g(i))) for i in range(self.order+1)]
    pdf = ROOT.RooAddPdf(f"laurent{self.order}", f"Laurent{self.order}", self.pdfs, params, True)
    pdf.fixCoefNormalization(x)
    return pdf
