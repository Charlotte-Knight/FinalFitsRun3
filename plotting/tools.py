import matplotlib.pyplot as plt
import mplhep
mplhep.style.use("CMS")
import os
import common.tools as tools
import logging
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

def histPlotTemplate(datahist, xlim):
  bin_centers, hist, uncert = tools.RooDataHist2Numpy(datahist, xlim=xlim)
  bin_width = bin_centers[1] - bin_centers[0]
  plt.errorbar(bin_centers, hist, xerr=bin_width/2, yerr=uncert, capsize=2, fmt='k.')

  plt.xlabel(r"$m_{\gamma\gamma}$ [GeV]")
  plt.ylabel(f"Events / ( {bin_width:.2g} GeV )")
  plt.ylim(bottom=0)
  cmslabel()
  return bin_width
