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
               order: int = 1, transforms = None, polys = None) -> None:
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
    self.init_param_bounds(bounds)
    self.init_transforms(transforms)
    self.init_polys(polys)
    self.init_params()
    self.init_x(x)
    self.init_roopdf()
    set_pre_postfix(self.roopdf, *self.params.values(), prefix=prefix, postfix=postfix)

  def __call__(self, xi):
    return utils.getVal(self.roopdf, self.x, xi)

  def get_final_shape_param_names(self):
    """Return all the names of the shape parameters used to initialize the pdf"""
    return [f"{name}{i+1}" for i in range(self.order) for name in self.default_bounds]
  
  def expand_config(self, config, default_config, require_match=False):
    """
    Helper method to expand the param bounds, transforms and polys config dictionaries.
    Is needed when for when wildcards are used.
    """
    config = config if config else default_config
    expanded_config = {}

    for name in self.get_final_shape_param_names():
      matches = [k for k in config if re.match(k, name)]
      if require_match:
        assert len(matches) == 1, f"No match or multiple matches for {name} in config"
      else:
        assert len(matches) <= 1, f"Multiple matches for {name} in config"
      
      expanded_config[name] = config[matches[0]] if matches else None
    return expanded_config
  
  def init_param_bounds(self, bounds: dict[str, tuple[float, float, float]]) -> None:
    """Initialize bounds for parameters of the pdf

    Args:
        bounds (dict[str, tuple[float, float, float]]): Dictionary with parameter names as keys and tuples with (default value, min, max) as values.
    """
    self.bounds = self.expand_config(bounds, self.default_bounds, require_match=True)

  def init_transforms(self, transforms):
    self.transforms = self.expand_config(transforms, self.default_transforms)
    for name, t in self.transforms.items():
      if t:
        t0 = t[0] if isinstance(t[0], ROOT.RooAbsReal) else ROOT.RooFit.RooConst(t[0])
        t1 = t[1] if isinstance(t[1], ROOT.RooAbsReal) else ROOT.RooFit.RooConst(t[1])
        self.transforms[name] = (t0, t1)

  def init_polys(self, polys):
    self.polys = self.expand_config(polys, {})
    for name, p in self.polys.items():
      if p and not isinstance(p[0], ROOT.RooAbsReal):
        self.polys[name][0] = ROOT.RooFit.RooConst(p[0])

  def init_params(self):
    """Initialize parameters (RooRealVars) of the pdf"""
    self.params = {}
    for name in self.get_final_shape_param_names():
      bounds = np.array(self.bounds[name])
      transform = self.transforms[name]
      poly = self.polys[name]

      if not (transform or poly):
        self.params[name] = ROOT.RooRealVar(name, name, *bounds)
        continue

      free_name = f"{name}_free" if transform else name
      bounds_transformed = (bounds - transform[0].getVal()) / transform[1].getVal() if transform else bounds

      if poly:
        poly_var, poly_order = poly

        pcn = lambda i: f"{free_name}_polycoeff{i}" #polycoeff name
        # # passing through the bounds during initialization not yet implemented
        polycoeff_bounds = [self.bounds[pcn(i)] if pcn(i) in self.bounds else (0.0, -0.02, 0.02) for i in range(poly_order+1)]
        polycoeff_bounds[1] = (bounds_transformed[0]/poly_var.getVal(), -0.02, 0.02)
        
        polycoeffs = [ROOT.RooRealVar(f"{free_name}_polycoeff{i}", f"{free_name}_polycoeff{i}", *polycoeff_bounds[i])
                      for i in range(poly_order+1)]
        
        formula_str = "+".join([f"@{i+1}*@0**{i}" for i in range(poly_order+1)])
        poly = ROOT.RooFormulaVar(free_name, free_name, formula_str, [poly_var]+polycoeffs)
        
        self.params[free_name] = poly
        for i, c in enumerate(polycoeffs):
          self.params[c.GetName()] = c
      else:
        self.params[free_name] = ROOT.RooRealVar(free_name, free_name, *bounds_transformed)
            
      if transform:
        self.params[name] = ROOT.RooFormulaVar(name, name, f"@0*@2+@1",
                                   [self.params[free_name], transform[0], transform[1]])
  
  def init_x(self, x) -> None:
    self.x = x
    if self.x_norm_factor != 1:
      self.x_norm = ROOT.RooFormulaVar("x_norm", "x_norm", f"@0*{self.x_norm_factor}", [x])
    else:
      self.x_norm = x

  def init_roopdf(self, name = None) -> None:
    """Initialize the RooAbsPdf object

    Args:
        x (ROOT.RooRealVar): the random variable (usually mass) of the pdf
    """
    name = name if name else f"{self.__class__.__name__}{self.order}"
    shape_params = [self.params[name] for name in self.get_final_shape_param_names()]
    self.roopdf = self.roopdf_constructor(name, name, self.x_norm, *shape_params)

  @property
  def poly_names(self):
    """Names of the params which are polynomials"""
    return [k for (k,v) in self.polys.items() if v]

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
    print("randomizing")
    rng = np.random.default_rng(seed)
    for p in self.free_params.values():
      print(p.GetName(), p.getMin(), p.getMax())
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
  # default_transforms = {
  #   "mean": [125, 1], 
  # }

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