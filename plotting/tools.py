import matplotlib.pyplot as plt
import mplhep
mplhep.style.use("CMS")
import os
import common.tools as tools
import logging
import numpy as np
log = logging.getLogger(__name__)

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
  bin_centers, hist, uncert = tools.RooDataHist2Numpy(datahist, xlim=xlim)
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
