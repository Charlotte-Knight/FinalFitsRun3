import logging

import numpy as np
import matplotlib.pyplot as plt
import mplhep
import ROOT

from finalfits import utils

log = logging.getLogger(__name__)
mplhep.style.use("CMS")
ROOT.gROOT.SetBatch(True)

def plotFit(datahist, x, pdf, params, savepath, xlim=None):
  log.info("Plotting fit")
  bin_centers, hist, uncert = utils.RooDataHist2Numpy(datahist, xlim=xlim)

  if xlim is None:
    xlim = (x.getMin(), x.getMax())

  bin_width = bin_centers[1] - bin_centers[0]
  plt.errorbar(bin_centers, hist, xerr=bin_width/2, yerr=uncert, capsize=2, fmt='k.')

  sf = datahist.sumEntries() * bin_width
  xi = np.linspace(xlim[0], xlim[1], 1000)
  plt.plot(xi, utils.getVal(pdf, x, xi)*sf)

  text = str(pdf.getTitle()) + " Fit"
  plt.text(0.05, 0.95, text, verticalalignment='top', transform=plt.gca().transAxes)
  text = ""
  for param in params:
    text += utils.textify(str(param.getTitle())) + r"$=%.2f \pm %.2f$"%(param.getVal(), param.getError()) + "\n"
  plt.text(0.06, 0.89, text, verticalalignment='top', transform=plt.gca().transAxes, fontsize='small')
  chi2 = ((hist-utils.getVal(pdf, x, bin_centers)*sf)**2 / uncert**2).sum() / len(hist) #chi2 per d.o.f
  plt.text(max(xi), max(hist+uncert), r"$\chi^2 / dof$=%.2f"%chi2, verticalalignment='top', horizontalalignment='right')
 
  plt.xlabel(r"$m_{\gamma\gamma}$")
  plt.ylim(bottom=0)
  utils.savefig(savepath)

def plotFitRoot(datahist, x, pdf, params, savepath, xlim=None):
  xc = x.Clone()
  if xlim is None:
    xlim = (x.getMin(), x.getMax())
  xc.setMin(xlim[0])
  xc.setMax(xlim[1])

  xframe = xc.frame()
  datahist.plotOn(xframe)
  pdf.plotOn(xframe)
  pdf.plotOn(xframe, Range="Full")
  c = ROOT.TCanvas("c", "c", 400, 400)
  xframe.Draw()
  c.SaveAs(f"{savepath}_root.pdf")

def plotFamily(datahist, x, results, savepath, blinded_regions=[], title=None):
  log.info("Plotting fTest fits")
  
  #xlim = (100, 180)
  xlim = (x.getMin(), x.getMax())
  bin_width = utils.histPlotTemplate(datahist, xlim, blinded_regions)

  sf = datahist.sumEntries() * bin_width
  xi = np.linspace(xlim[0], xlim[1], 1000)
  for res in results:
    label = f"{title} {res['dof']}: " + r"$p_{ftest}=%.2f$, "%res["ftest_pval"] + r"$p_{gof}=%.2f$ "%res["gof_pval"]
    plt.plot(xi, utils.getVal(res["pdf"].roopdf, x, xi)*sf, label=label)

  plt.legend()
  utils.savefig(savepath)
  
def plotEnvelope(datahist, x, results, savepath, blinded_regions=[]):
  log.info("Plotting envelope")

  #xlim = (100, 180)
  xlim = (x.getMin(), x.getMax())
  bin_width = utils.histPlotTemplate(datahist, xlim, blinded_regions)
  
  sf = datahist.sumEntries() * bin_width
  xi = np.linspace(xlim[0], xlim[1], 1000)

  gofs = [res["gof_pval"] for family in results.keys() for res in results[family]]
  best_gof_index = gofs.index(max(gofs))

  for family in results.keys():
    for res in results[family]:
      label = f"{family} {res['dof']}"
      plt.plot(xi, utils.getVal(res["pdf"].roopdf, x, xi)*sf, label=label)

  legend = plt.legend()
  handle = legend.get_texts()[best_gof_index]
  handle.set_weight("bold")
  text = handle.get_text() + "\n(Best fit)"
  handle.set_text(text)

  utils.savefig(savepath)