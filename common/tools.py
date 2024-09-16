import ROOT
import numpy as np
import logging
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
  dataset = w.data("data")
  return x, dataset

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

def checkFunctionAtBounds(f):
  for var in f.vars:
    if np.isclose(var.getVal(), var.getMin(), rtol=0.01) or np.isclose(var.getVal(), var.getMax(), rtol=0.01):
      log.warning(f"Parameter {var.GetName()} from function {f.pdf.GetName()} is at its bounds")
      log.warning(f"{var.GetName()}={var.getVal()}, low={var.getMin()}, high={var.getMax()}")