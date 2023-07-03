import numpy as np
import matplotlib.pyplot as plt
import mplhep
mplhep.style.use("CMS")

import plotting.tools as ptools

title_dict = {
      #"mean":r"$\bar{m}_{\gamma\gamma}$",
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

def plotFit(datahist, x, pdf, params, savepath, xlim=None):
  bin_centers, hist, uncert = ptools.RooDataHist2Numpy(datahist, xlim=xlim)

  if xlim is None:
    xlim = (x.getMin(), x.getMax())

  bin_width = bin_centers[1] - bin_centers[0]
  plt.errorbar(bin_centers, hist, xerr=bin_width/2, yerr=uncert, capsize=2, fmt='k.')

  sf = datahist.sumEntries() * bin_width
  xi = np.linspace(xlim[0], xlim[1], 1000)
  plt.plot(xi, ptools.getVal(pdf, x, xi)*sf)

  text = str(pdf.getTitle()) + " Fit"
  plt.text(0.05, 0.95, text, verticalalignment='top', transform=plt.gca().transAxes)
  text = ""
  for param in params:
    text += textify(str(param.getTitle())) + r"$=%.2f \pm %.2f$"%(param.getVal(), param.getError()) + "\n"
  plt.text(0.06, 0.89, text, verticalalignment='top', transform=plt.gca().transAxes, fontsize='small')
  chi2 = ((hist-ptools.getVal(pdf, x, bin_centers)*sf)**2 / uncert**2).sum() / len(hist) #chi2 per d.o.f
  plt.text(max(xi), max(hist+uncert), r"$\chi^2 / dof$=%.2f"%chi2, verticalalignment='top', horizontalalignment='right')
 
  plt.xlabel(r"$m_{\gamma\gamma}$")
  plt.ylim(bottom=0)
  plt.savefig(f"{savepath}.png")
  plt.savefig(f"{savepath}.pdf")

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