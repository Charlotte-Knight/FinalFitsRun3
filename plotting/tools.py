import ROOT
import numpy as np
import matplotlib.pyplot as plt
import mplhep
mplhep.style.use("CMS")

def RooDataHist2Numpy(datahist, xlim=None):
  nBins = datahist.get()[0].getBins()
  
  bin_centers = []
  hist = []

  for i in range(nBins):
    bin_center = datahist.get(i)[0].getVal()
    if (xlim is None) or (xlim[0] <= bin_center <= xlim[1]):
      bin_centers.append(datahist.get(i)[0].getVal())
      hist.append(datahist.weight(i))

  bin_centers = np.array(bin_centers)
  hist = np.array(hist)
  uncert = np.sqrt(hist)
  #uncert[uncert==0] = uncert[uncert!=0].min()
  uncert[uncert==0] = 1

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