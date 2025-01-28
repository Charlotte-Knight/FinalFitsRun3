import os
import logging

import numpy as np
import matplotlib.pyplot as plt
import ROOT
import mplhep

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
      uncert.append(datahist.weightError(etype=ROOT.RooAbsData.ErrorType.SumW2))

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

def getNBinsFitted(x, fit_ranges):
  bin_boundaries = np.linspace(x.getMin(), x.getMax(), x.getBins()+1)
  bin_centers = (bin_boundaries[:-1] + bin_boundaries[1:]) / 2
  
  inside_ranges = np.zeros_like(bin_centers, dtype=bool)
  for r in fit_ranges:
    inside_ranges = inside_ranges | ((bin_centers > r[0]) & (bin_centers < r[1]))

  nbins_fitted = sum(inside_ranges)
  return nbins_fitted


title_dict = {
      "mean":r"$\mu$",
      "sigma":r"$\sigma$",
      "sigmaLR":r"$\sigma$",
      "alphaL":r"$\alpha_L$",
      "nL":r"$n_L$",
      "alphaR":r"$\alpha_R$",
      "nR":r"$n_R$",
      "c":r"$c$"
      }

def textify(title):
  if (title not in title_dict) and (title[:-1] not in title_dict):
    return title
  else:
    if title[-1].isdigit():
      idx = int(title[-1])
      title = title[:-1]
    else:
      idx = None

    title = title_dict[title]
    if idx is not None:
      title = title[:-1] + "_{%d}$"%idx
    
    return title
  
def savefig(savepath, extensions=["png", "pdf"], keep=False):
  directory = "/".join(savepath.split("/")[:-1])
  os.makedirs(directory, exist_ok=True)
  for extension in extensions:
    log.info(f"Saving figure to {savepath}.{extension}")
    plt.savefig(f"{savepath}.{extension}")
  if not keep:
    plt.clf()

def cmslabel():
  mplhep.cms.label("Work in Progress", data=True, lumi=138, com=13.6)

def histPlotTemplate(datahist, xlim, blinded_regions=[]):
  bin_centers, hist, uncert = RooDataHist2Numpy(datahist, xlim=xlim)
  bin_width = bin_centers[1] - bin_centers[0]

  blinded_idx = []
  for region in blinded_regions:
    for idx, bin_center in enumerate(bin_centers):
      if float(region.split(",")[0]) < bin_center < float(region.split(",")[1]):
        blinded_idx.append(idx)
  
  s = ~np.isin(np.arange(len(bin_centers)), blinded_idx)

  plt.errorbar(bin_centers[s], hist[s], xerr=bin_width/2, yerr=uncert[s], capsize=2, fmt='k.')

  plt.xlabel(r"$m_{\gamma\gamma}$ [GeV]")
  plt.ylabel(f"Events / ( {bin_width:.2g} GeV )")
  plt.ylim(bottom=0)
  cmslabel()
  return bin_width

finalfits_verbose_dict = {
  -2: "CRITICAL",
  -1: "ERROR",
  0:  "WARNING",
  1:  "INFO",
  2:  "DEBUG"
}

roofit_verbose_dict = {
  -2: "FATAL",
  -1: "ERROR",
  0:  "WARNING",
  1:  "PROGRESS",
  2:  "INFO",
  3:  "DEBUG"
}

def comma_separated_two_tuple(string):
  numbers = string.split(",")
  if len(numbers) != 2:
    raise TypeError("Must be two numbers")
  return tuple(map(float, numbers))

def addLoggingArguments(parser):
  parser.add_argument("--verbose", "-v", type=int, default=1, choices=range(-2,3),
                      help="Set verbosity level for finalFits scripts: %s"%(", ".join(f"{key}={value}" for key,value in finalfits_verbose_dict.items())))
  parser.add_argument("--roofit-verbose", type=int, default=0, choices=range(-2,4),
                      help="Set verbosity level for RooFit: %s"%(", ".join(f"{key}={value}" for key,value in roofit_verbose_dict.items())))

def applyLoggingArguments(args):
  ROOT.RooMsgService.instance().setGlobalKillBelow(getattr(ROOT.RooFit, roofit_verbose_dict[args.roofit_verbose]))
  logging.basicConfig(level=getattr(logging, finalfits_verbose_dict[args.verbose]), format=('%(name)-20s: %(levelname)-8s %(message)s'))