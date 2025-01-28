"""
Module that helps build a signal model that can be simulatenously fit to different mass points.

Some of the model parameters will be built as polynomials which are functions of the mass.
"""

import logging
import re
from typing import Union

import ROOT

from finalfits import pdfs

log = logging.getLogger(__name__)

def get_roocategory(category_names):
  cat = ROOT.RooCategory("cat", "cat")
  for cat_name in category_names:
    cat.defineType(cat_name)
  return cat

def get_sim_model(x: ROOT.RooRealVar, MH: ROOT.RooRealVar, mass_points: list[int], pdf_name: str, order: int):
  pdf = getattr(pdfs, pdf_name)(x, order=order, 
                                transforms={"mean*": [MH, 1]},
                                polys={"mean*|sigma*": [MH, 1]})
  
  mass_pdfs = [getattr(pdfs, pdf_name)(x, order=order, postfix=f"m{m}",
                                       transforms={"mean*": [m, 1]}, polys={"mean*|sigma*": [m, 1]})
               for m in mass_points]
  
  # want to share all parameters except the polynomials
  shared_param_names = self.params.keys() - self.poly_names
  shared_params = {name: self.params[name] for name in shared_param_names}
  for mpdf in mass_pdfs:
    mpdf.overide_params(shared_params)