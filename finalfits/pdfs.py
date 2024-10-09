"""
Module that contains the pdfs (probability density functions) that can be used by finalfits.

Every pdf is a class that inherits from FinalFitsPdf which wraps around RooFit functionality.
A FinalFitsPdf instance contains the RooFit pdf object and its parameters. 
"""

import logging
from typing import Optional

import numpy as np
import ROOT

log = logging.getLogger(__name__)
available_pdfs = ["Gaussian", "DCB", "Bernstein", "Exponential", "Power", "ExpPoly", "Laurent"]

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
  roopdf_constructor = ROOT.RooAbsPdf # to be overwritten by subclass
  default_bounds = {"param1": [1, 0, 2], } #... to be overwritten by subclass
  x_norm_factor = 1 # factor to multiply x by (can be useful to get x~1)
  max_order = 5
    
  def __init__(self, x: ROOT.RooRealVar, prefix: str = "", postfix: str = "",
               bounds: Optional[dict[str, tuple[float, float, float]]] = None,
               order: int = 1) -> None:
    """
    Args:
        x (ROOT.RooRealVar): the random variable (usually mass) of the pdf
        prefix (str, optional): string to append to start of all RooFit object(s) name/title. Defaults to "".
        postfix (str, optional): string to append to start of all RooFit object(s) name/title. Defaults to "".
        bounds (Optional[dict[str, tuple[float, float, float]]], optional): Parameter bounds. Dictionary with parameter names as keys and tuples with (default value, min, max) as values. Defaults to None.
        order (int, optional): Order of the pdf. Defaults to 1.
    """
    if order > self.max_order:
      raise ValueError(f"Order of {self.__class__.__name__} is too high. Max order is {self.max_order}.")
    self.order = order
    self.init_parameter_bounds(bounds)
    self.init_params()
    self.init_x(x)
    self.init_roopdf()
    set_pre_postfix(self.roopdf, *self.params, prefix=prefix, postfix=postfix)

  def __str__(self) -> str:
    s = f"{self.__class__.__name__} of order {self.order} with parameters:\n"
    for p in self.params:
      s += f"{p.GetName()} = {p.getVal():.2f}, "
    return s

  def init_parameter_bounds(self, bounds: dict[str, tuple[float, float, float]]) -> None:
    """Initialize bounds for parameters of the pdf

    Args:
        bounds (dict[str, tuple[float, float, float]]): Dictionary with parameter names as keys and tuples with (default value, min, max) as values.
    """
    self.bounds = bounds if bounds is not None else self.default_bounds

  def init_params(self):
    """Initialize parameters (RooRealVars) of the pdf"""
    if self.order == 1:
      self.params = [ROOT.RooRealVar(name, name, *self.bounds[name]) for name in self.bounds]
    else:
      self.params = [ROOT.RooRealVar(f"{name}{i+1}", f"{name}{i+1}", *self.bounds[name]) 
                     for i in range(self.order) for name in self.bounds]

  def init_x(self, x) -> None:
    self.x = x
    if self.x_norm_factor != 1:
      self.x_norm = ROOT.RooFormulaVar("x_norm", "x_norm", f"@0*{self.x_norm_factor}", [x])
    else:
      self.x_norm = x

  def init_roopdf(self) -> None:
    """Initialize the RooAbsPdf object

    Args:
        x (ROOT.RooRealVar): the random variable (usually mass) of the pdf
    """
    name = f"{self.__class__.__name__}{self.order}"
    self.roopdf = self.roopdf_constructor(name, name, self.x_norm, *self.params)

  def randomize_params(self, skip_constant: bool = True, seed: bool = None) -> None:
    """randomize the parameters of the pdf

    Args:
        skip_constant (bool, optional): skip variables which are set to be constant in RooFit. Defaults to True.
        seed (bool, optional): random seed. Defaults to None.
    """
    np.random.seed(seed)
    for v in self.params:
      if (not v.isConstant()) or (not skip_constant):
        v.setVal(np.random.uniform(v.getMin(), v.getMax()))

  def get_dof(self) -> int:
    """Get the degrees of freedom of the pdf. This is the number of parameters of the pdf.

    Returns:
        int: degrees of freedom of the pdf
    """
    return len(self.params)

  def check_bounds(self) -> None:
    """Check if any of the parameters are at their bounds and log a warning if they are."""
    for p in self.params:
      if np.isclose(p.getVal(), p.getMin(), rtol=0.01) or np.isclose(p.getVal(), p.getMax(), rtol=0.01):
        log.warning("Parameter %s from pdf %s is at its bounds", p.GetName(), self.roopdf.GetName())
        log.warning("%s=%s, low=%f, high=%f", p.GetName(), p.getVal(), p.getMin(), p.getMax())

class FinalFitsPdfSum(FinalFitsPdf):
  """Base class for pdfs that are sums of other pdfs."""
  def init_params(self) -> None:
    super().init_params()
    self.coeffs = [ROOT.RooRealVar(f"c{i}", f"c{i}", 1/self.order, 0.0, 1.0) for i in range(self.order-1)]
    self.params += self.coeffs    
    
  def init_roopdf(self) -> None:
    if self.order == 1:
      super().init_roopdf()
    else:
      self.pdfs = []
      name = f"{self.__class__.__name__}{self.order}" 
      
      for i in range(self.order):
        subname = f"{name}component{i+1}"
        params_subset = self.params[i*len(self.bounds):(i+1)*len(self.bounds)]
        self.pdfs.append(self.roopdf_constructor(subname, subname, self.x_norm, *params_subset))
        
      self.roopdf = ROOT.RooAddPdf(name, name, self.pdfs, self.coeffs, True)
      self.roopdf.fixCoefNormalization(self.x_norm)

class Gaussian(FinalFitsPdfSum):
  """Sum of Gaussians pdf"""
  roopdf_constructor = ROOT.RooGaussian
  default_bounds = {
    "mean": [125, 123, 127],
    "sigma": [1.5, 1, 5]
  }

class DCB(FinalFitsPdf):
  """Double Crystal Ball pdf"""
  roopdf_constructor = ROOT.RooCrystalBall
  default_bounds = {
    "mean": [125, 120, 130],
    "sigmaLR": [1.3, 1, 5],
    "alphaL": [1, 0.1, 5],
    "nL": [5, 0.1, 20],
    "alphaR": [1, 0.1, 5],
    "nR": [20, 0.1, 20]
  }
  max_order = 1

class Exponential(FinalFitsPdfSum):
  """Sum of exponentials pdf"""
  roopdf_constructor = ROOT.RooExponential
  default_bounds = {"a": [-0.05, -0.5, 0]}

class Power(FinalFitsPdfSum):  
  """Sum of power laws pdf"""
  roopdf_constructor = ROOT.RooPower
  default_bounds = {"a": [-1, -5, 0]}

class Bernstein(FinalFitsPdf):
  """Bernstein polynomial pdf"""
  def roopdf_constructor(self, name: str, title: str, x: ROOT.RooRealVar, *params: ROOT.RooRealVar) -> ROOT.RooAbsPdf:
    return ROOT.RooBernstein(name, title, x, [ROOT.RooFit.RooConst(1.0), *params])
  default_bounds = {"a": [0.1, 0, 10]}

class ExpPoly(FinalFitsPdf):
  """Exponential polynomial pdf"""
  def roopdf_constructor(self, name: str, title: str, x: ROOT.RooRealVar, *params: ROOT.RooRealVar) -> ROOT.RooAbsPdf:
    return ROOT.RooExpPoly(name, title, x, ROOT.RooArgList(*params))
  default_bounds = {"a": [-0.5, -1, 0]}
  x_norm_factor = 0.01

class Laurent(FinalFitsPdf):
  """Laurent polynomial pdf"""
  default_bounds = {"a": [0.5, 0, 1]}
  
  def init_roopdf(self) -> None:
    def g(i):
      return sum([(-1)**j * j for j in range(i+1)])
    
    name = f"{self.__class__.__name__}{self.order}"
    self.pdfs = [ROOT.RooPower(f"{name}component{i}", f"{name}component{i}", self.x_norm, 
                               ROOT.RooFit.RooConst(-4+g(i))) for i in range(self.order+1)]
    self.roopdf = ROOT.RooAddPdf(name, name, self.pdfs, self.params, True)
    self.roopdf.fixCoefNormalization(self.x_norm)