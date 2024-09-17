import ROOT
import numpy as np
import logging
import common.functions as functions
log = logging.getLogger(__name__)

def RooDataHist2Numpy(datahist, xlim=None):
  nBins = datahist.get()[0].getBins()
  
  bin_centers = []
  hist = []
  uncert = []

  for i in range(nBins):
    bin_center = datahist.get(i)[0].getVal()
    if (xlim is None) or (xlim[0] <= bin_center <= xlim[1]):
      bin_centers.append(datahist.get(i)[0].getVal())
      hist.append(datahist.weight(i))
      uncert.append(datahist.weightError())

  bin_centers = np.array(bin_centers)
  hist = np.array(hist)
  uncert = np.array(uncert)

  return bin_centers, hist, uncert

def getVal(pdf, xvar, xval):
  if hasattr(xval, "__len__"):
    vals = []
    for xi in xval:
      xvar.setVal(xi)
      vals.append(pdf.getVal())
    val = np.array(vals)
  else:
    xvar.setVal(xval)
    val = pdf.getVal()
  return val / pdf.createIntegral(xvar).getVal()

def readEvents(filename):
  log.info(f"Loading workspace from {filename}")
  f = ROOT.TFile(filename)
  w = f.Get("w")
  x = w.var("x")
  data = w.data("data")
  return x, data

def setRanges(x, fit_ranges):
  x.setRange("Full", x.getMin(), x.getMax())
  for i, r in enumerate(fit_ranges):
    low, high = r.split(",")
    low, high = float(low), float(high)
    x.setRange(f"range{i}", low, high)
  fit_ranges_str = ",".join([f"range{i}" for i, r in enumerate(fit_ranges)])
  return fit_ranges_str

def getNBinsFitted(x, fit_ranges):
  bin_boundaries = np.linspace(x.getMin(), x.getMax(), x.getBins()+1)
  bin_centers = (bin_boundaries[:-1] + bin_boundaries[1:]) / 2
  inside_ranges = np.zeros_like(bin_centers, dtype=bool)

  for i, r in enumerate(fit_ranges):
    low, high = r.split(",")
    low, high = float(low), float(high)
    inside_ranges = inside_ranges | ((bin_centers > low) & (bin_centers < high))
  nbins_fitted = sum(inside_ranges)
  return nbins_fitted

def checkPdfAtBounds(pdf, datahist):
  vars = pdf.getParameters(datahist)
  for var in vars:
    if np.isclose(var.getVal(), var.getMin(), rtol=0.01) or np.isclose(var.getVal(), var.getMax(), rtol=0.01):
      log.warning(f"Parameter {var.GetName()} from function {pdf.GetName()} is at its bounds")
      log.warning(f"{var.GetName()}={var.getVal()}, low={var.getMin()}, high={var.getMax()}")

def robustFit(pdf, datahist, fit_ranges_str, n_fits=4, recursive=True):
  NLLs = []
  vars_vals = []  

  log.info(f"Fitting {pdf.GetName()} to {datahist.GetName()}. Doing {n_fits} fits from random initialisations.")
  for i in range(n_fits):
    functions.randomiseVars(pdf.getParameters(datahist))
    r = pdf.fitTo(datahist, Range=fit_ranges_str, PrintLevel=-1, Save=True)
    NLLs.append(r.minNll())
    vars_vals.append([var.getVal() for var in pdf.getParameters(datahist)])

  best_vars_vals = vars_vals[np.argmin(NLLs)]
  for i, val in enumerate(best_vars_vals):
    pdf.getParameters(datahist)[i].setVal(val)

  max_diff = 0.01
  max_n_fits = 1024

  NLL1 = min(NLLs[:n_fits//2])
  NLL2 = min(NLLs[n_fits//2:])
  NLL_diff = abs(NLL1 - NLL2)
  log.debug(f"Difference in minimum NLL from first and second half of fits is {NLL_diff}")

  if NLL_diff > max_diff:
    log.warning(f"Fit is unstable when starting from {n_fits} random initialisations.")
    
    if recursive and n_fits < max_n_fits:
      n_fits_more = n_fits * 2
      robustFit(pdf, datahist, fit_ranges_str, n_fits_more)
    else:
      raise Exception(f"Fit is unstable. Tried {n_fits} from random places and did not find acceptable convergence. Halted because we reached maximum number of fits ({max_n_fits}).")

  checkPdfAtBounds(pdf, datahist)