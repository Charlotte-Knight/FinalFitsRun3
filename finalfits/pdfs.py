"""
Module that contains the pdfs (probability density functions) that can be used by finalfits.

Every pdf is a class that inherits from FinalFitsPdf which wraps around RooFit functionality.
A FinalFitsPdf instance contains the RooFit pdf object and its parameters. 
"""

import logging
from typing import Optional
import re

import numpy as np
import ROOT

from finalfits import utils

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
  default_bounds = {} # {"param1": [1, 0, 2], } #... to be overwritten by subclass
  default_transforms = {} # {"param1": [0, 1], } #... to be overwritten by subclass
  x_norm_factor = 1 # factor to multiply x by (can be useful to get x~1)
  max_order = 5
    
  def __init__(self, x: ROOT.RooRealVar, prefix: str = "", postfix: str = "",
               bounds: Optional[dict[str, tuple[float, float, float]]] = None,
               order: int = 1, transforms = None) -> None:
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
    self.init_transforms(transforms)
    self.init_parameter_bounds(bounds)
    self.init_params()
    self.init_x(x)
    self.init_roopdf()
    set_pre_postfix(self.roopdf, *self.params.values(), prefix=prefix, postfix=postfix)

  def __call__(self, xi):
    return utils.getVal(self.roopdf, self.x, xi)

  # def __str__(self) -> str:
  #   s = f"{self.__class__.__name__} of order {self.order} with parameters:\n"
  #   for p in self.params:
  #     s += f"{p.GetName()} = {p.getVal():.2f}, "
  #   return s

  def init_transforms(self, transforms):
    transforms = transforms if transforms is not None else self.default_transforms
    self.transforms = {}
    for name, t in transforms.items():
      t0 = t[0] if isinstance(t[0], ROOT.RooAbsReal) else ROOT.RooFit.RooConst(t[0])
      t1 = t[1] if isinstance(t[1], ROOT.RooAbsReal) else ROOT.RooFit.RooConst(t[1])
      self.transforms[name] = (t0, t1)
    
  def init_parameter_bounds(self, bounds: dict[str, tuple[float, float, float]]) -> None:
    """Initialize bounds for parameters of the pdf

    Args:
        bounds (dict[str, tuple[float, float, float]]): Dictionary with parameter names as keys and tuples with (default value, min, max) as values.
    """
    self.bounds = bounds if bounds is not None else self.default_bounds    

  def get_transform(self, name):
    matches = []
    for k, v in self.transforms.items():
      if re.match(k, name):
        matches.append(v)
    assert len(matches) <= 1, f"Multiple matches for {name} in {self.transforms}"
    return matches[0] if matches else None

  def init_param(self, name, i=None):
    var_name = name if i is None else f"{name}{i+1}"
    bounds = np.array(self.bounds[name])
    
    t = self.get_transform(name)
    if t:
      # useful to keep track later
      self.transforms[name] = t
      self.transforms[var_name] = t
      
      bounds_transformed = sorted((bounds - t[0].getVal()) / t[1].getVal())
      
      var_free_name = f"{var_name}_free"
      self.params[var_free_name] = ROOT.RooRealVar(var_free_name, var_free_name, *bounds_transformed)
      self.params[var_name] = ROOT.RooFormulaVar(var_name, var_name, f"@0*@2+@1", 
                                               [self.params[var_free_name], t[0], t[1]])
    else:
      self.params[var_name] = ROOT.RooRealVar(var_name, var_name, *bounds)

  def init_params(self):
    """Initialize parameters (RooRealVars) of the pdf"""
    self.params = {}
    if self.order == 1:
      for name in self.bounds:
        self.init_param(name)
    else:
      for i in range(self.order):
        for name in self.bounds:
          self.init_param(name, i)
  
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

    if self.order == 1:
      shape_params = [self.params[name] for name in self.default_bounds]
    else:
      shape_params = [self.params[f"{name}{i+1}"] for i in range(self.order) for name in self.default_bounds]

    self.roopdf = self.roopdf_constructor(name, name, self.x_norm, *shape_params)

  def overide_params(self, params: dict[str, ROOT.RooAbsReal]) -> None:
    for name, p in params.items():
      self.params[name] = p    

  @property
  def free_params(self):
    return {k: v for (k, v) in self.params.items()
            if isinstance(v, ROOT.RooRealVar) and not v.isConstant()}

  @property
  def free_params_vals(self):
    return {k: v.getVal() for (k, v) in self.free_params.items()}
  
  @free_params_vals.setter
  def free_params_vals(self, vals):
    for k, v in vals.items():
      self.params[k].setVal(v)
      
  @property
  def final_params(self):
    return {k: v for (k, v) in self.params.items() if "_free" not in k}
  
  @property
  def final_params_vals(self):
    return {k: v.getVal() for (k, v) in self.final_params.items()}
  
  @property
  def final_params_errs(self):
    final_params_errs = {}
    for name, p in self.final_params.items():
      if isinstance(p, ROOT.RooRealVar):
        err = p.getError()
      else:
        err = self.params[name+"_free"].getError() * self.transforms[name][1].getVal()

      final_params_errs[name] = err
    return final_params_errs
      
  def randomize_params(self, seed: int = None) -> None:
    """randomize the parameters of the pdf

    Args:
        seed (bool, optional): random seed. Defaults to None.
    """
    rng = np.random.default_rng(seed)
    for p in self.free_params.values():
      p.setVal(rng.uniform(p.getMin(), p.getMax()))

  def get_dof(self) -> int:
    """Get the degrees of freedom of the pdf. This is the number of free parameters of the pdf.

    Returns:
        int: degrees of freedom of the pdf
    """
    return len(self.free_params)

  def check_bounds(self) -> None:
    """Check if any of the parameters are at their bounds and log a warning if they are."""
    for p in self.free_params.values():
      if np.isclose(p.getVal(), p.getMin(), rtol=0.01) or np.isclose(p.getVal(), p.getMax(), rtol=0.01):
        log.warning("Parameter %s from pdf %s is at its bounds", p.GetName(), self.roopdf.GetName())
        log.warning("%s=%s, low=%f, high=%f", p.GetName(), p.getVal(), p.getMin(), p.getMax())
        
class FinalFitsPdfSum(FinalFitsPdf):
  """Base class for pdfs that are sums of other pdfs."""
  def init_params(self) -> None:
    super().init_params()
    self.coeffs = [ROOT.RooRealVar(f"c{i}", f"c{i}", 1/self.order, 0.0, 1.0) 
                   for i in range(self.order-1)]
    self.params.update({c.GetName(): c for c in self.coeffs})
    
  def init_roopdf(self) -> None:
    if self.order == 1:
      super().init_roopdf()
    else:
      self.pdfs = []
      name = f"{self.__class__.__name__}{self.order}"
      
      for i in range(self.order):
        subname = f"{name}component{i+1}"
        print(self.params)
        params_subset = [self.params[f"{param_name}{i+1}"] for param_name in self.default_bounds]
        self.pdfs.append(self.roopdf_constructor(subname, subname, self.x_norm, *params_subset))
        
      self.roopdf = ROOT.RooAddPdf(name, name, self.pdfs, self.coeffs, True)
      self.roopdf.fixCoefNormalization(self.x_norm)

class Gaussian(FinalFitsPdfSum):
  roopdf_constructor = ROOT.RooGaussian
  """Sum of Gaussians pdf"""
  default_bounds = {
    "mean": [125, 120, 130],
    "sigma": [1.5, 1, 5]
  }
  default_transforms = {
    "mean": [125, 1], 
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
    "nR": [10, 0.1, 20]
  }
  default_transforms = {
    "mean": [125, 1], 
  }
  max_order = 1

class Exponential(FinalFitsPdfSum):
  """Sum of exponentials pdf"""
  roopdf_constructor = ROOT.RooExponential
  default_bounds = {"a": [-0.02, -0.05, 0]}

class Power(FinalFitsPdfSum):  
  """Sum of power laws pdf"""
  roopdf_constructor = ROOT.RooPower
  default_bounds = {"a": [-1, -2, 0]}

class Bernstein(FinalFitsPdf):
  """Bernstein polynomial pdf"""
  def roopdf_constructor(self, name: str, title: str, x: ROOT.RooRealVar, *params: ROOT.RooRealVar) -> ROOT.RooAbsPdf:
    return ROOT.RooBernstein(name, title, x, [ROOT.RooFit.RooConst(1.0), *params])
  default_bounds = {"a": [0.1, 0, 1]}

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
    self.roopdf = ROOT.RooAddPdf(name, name, self.pdfs, list(self.params.values()), True)
    self.roopdf.fixCoefNormalization(self.x_norm)