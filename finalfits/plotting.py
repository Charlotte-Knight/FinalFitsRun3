import logging

import numpy as np
import matplotlib.pyplot as plt
import mplhep
import ROOT

from finalfits import utils, pdfs

log = logging.getLogger(__name__)
mplhep.style.use("CMS")
ROOT.gROOT.SetBatch(True)

def plotHist(datahist, savepath=None, xlim=None):
  log.info("Plotting histogram")
  bin_centers, hist, uncert = utils.RooDataHist2Numpy(datahist, xlim=xlim)

  bin_width = bin_centers[1] - bin_centers[0]
  plt.errorbar(bin_centers, hist, xerr=bin_width/2, yerr=uncert, capsize=2, fmt='k.')
  plt.xlabel(r"$m_{\gamma\gamma}$")
  plt.ylim(bottom=0)
  if savepath:
    utils.savefig(savepath)

  return bin_centers, hist, uncert
  
def plotFit(datahist, pdf, x, savepath=None, xlim=None):
  if isinstance(pdf, pdfs.FinalFitsPdf):
    roopdf = pdf.roopdf
  elif isinstance(pdf, ROOT.RooAbsPdf):
    roopdf = pdf
  else:
    raise ValueError("pdf must be either a FinalFitsPdf or a RooAbsPdf")
  
  log.info("Plotting fit")
  if xlim is None:
    xlim = (x.getMin(), x.getMax())
  
  bin_centers, hist, uncert = plotHist(datahist, None, xlim)
  bin_width = bin_centers[1] - bin_centers[0]
  sf = datahist.sumEntries() * bin_width
  xi = np.linspace(xlim[0], xlim[1], 1000)
  plt.plot(xi, utils.getVal(roopdf, x, xi)*sf)

  text = str(roopdf.getTitle()) + " Fit"
  plt.text(0.05, 0.95, text, verticalalignment='top', transform=plt.gca().transAxes)
  chi2 = ((hist-utils.getVal(roopdf, x, bin_centers)*sf)**2 / uncert**2).sum() / len(hist) #chi2 per d.o.f
  plt.text(max(xi), max(hist+uncert), r"$\chi^2 / dof$=%.2f"%chi2, verticalalignment='top', horizontalalignment='right')
  
  if isinstance(pdf, pdfs.FinalFitsPdf):
    vals = pdf.final_params_vals
    errs = pdf.final_params_errs
    text = ""
    for name in vals:
      text += utils.textify(name) + r"$=%.2f \pm %.2f$"%(vals[name], errs[name]) + "\n"
    plt.text(0.06, 0.89, text, verticalalignment='top', transform=plt.gca().transAxes, fontsize='small')
  
  if savepath:
    utils.savefig(savepath)

def plotFitRoot(datahist, pdf, savepath, xlim=None):
  xc = pdf.x.Clone()
  if xlim is None:
    xlim = (pdf.x.getMin(), pdf.x.getMax())
  xc.setMin(xlim[0])
  xc.setMax(xlim[1])

  xframe = xc.frame()
  datahist.plotOn(xframe)
  pdf.roopdf.plotOn(xframe)
  pdf.roopdf.plotOn(xframe, Range="Full")
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