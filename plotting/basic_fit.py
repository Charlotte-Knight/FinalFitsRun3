import numpy as np
import matplotlib.pyplot as plt
import mplhep
mplhep.style.use("CMS")
import logging
log = logging.getLogger(__name__)

import common.tools as tools
import plotting.tools as ptools

import ROOT
ROOT.gROOT.SetBatch(True)

def plotFit(datahist, x, pdf, params, savepath, xlim=None):
  log.info("Plotting fit")
  bin_centers, hist, uncert = tools.RooDataHist2Numpy(datahist, xlim=xlim)

  if xlim is None:
    xlim = (x.getMin(), x.getMax())

  bin_width = bin_centers[1] - bin_centers[0]
  plt.errorbar(bin_centers, hist, xerr=bin_width/2, yerr=uncert, capsize=2, fmt='k.')

  sf = datahist.sumEntries() * bin_width
  xi = np.linspace(xlim[0], xlim[1], 1000)
  plt.plot(xi, tools.getVal(pdf, x, xi)*sf)

  text = str(pdf.getTitle()) + " Fit"
  plt.text(0.05, 0.95, text, verticalalignment='top', transform=plt.gca().transAxes)
  text = ""
  for param in params:
    text += ptools.textify(str(param.getTitle())) + r"$=%.2f \pm %.2f$"%(param.getVal(), param.getError()) + "\n"
  plt.text(0.06, 0.89, text, verticalalignment='top', transform=plt.gca().transAxes, fontsize='small')
  chi2 = ((hist-tools.getVal(pdf, x, bin_centers)*sf)**2 / uncert**2).sum() / len(hist) #chi2 per d.o.f
  plt.text(max(xi), max(hist+uncert), r"$\chi^2 / dof$=%.2f"%chi2, verticalalignment='top', horizontalalignment='right')
 
  plt.xlabel(r"$m_{\gamma\gamma}$")
  plt.ylim(bottom=0)
  ptools.savefig(savepath)

def plotFitRoot(datahist, x, pdf, params, savepath, xlim=None):
  xframe = x.frame()
  datahist.plotOn(xframe)
  pdf.plotOn(xframe, ROOT.RooFit.Range("Full"))
  c = ROOT.TCanvas("c", "c", 400, 400)
  xframe.Draw()
  c.SaveAs(f"{savepath}_root.pdf")

if __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser(
                    prog='Signal Fit Plotter',
                    description='Plots signal fits')
  parser.add_argument("in_file")
  parser.add_argument("out_file")
  parser.add_argument("--xlim", nargs=2, type=float, default=(115, 135))
  args = parser.parse_args()

  import ROOT
  f = ROOT.TFile(args.in_file)
  wsig = f.Get("wsig")

  x = wsig.var("x")
  datahist = wsig.data("datahist")
  pdf = wsig.pdf("gauss_cat0")
  params = pdf.getParameters(x)

  plotFit(datahist, x, pdf, params, args.out_file, xlim=args.xlim)