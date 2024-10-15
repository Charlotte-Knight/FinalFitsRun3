"""
Module containing tools to perform fits (signal or background) to datasets (data or MC).
This includes:
- Generic robust fitting function
- Ftest (including goodness-of-fit tests)
"""
import logging

import numpy as np
import ROOT
from finalfits import pdfs
from finalfits import utils
from finalfits import plotting

log = logging.getLogger(__name__)

def prepare_ranges(x, fit_ranges):
  if fit_ranges == ():
    fit_ranges = ((x.getMin(), x.getMax()), )

  x.setRange("Full", x.getMin(), x.getMax())

  fit_ranges_dict = {f"range{i}": r for i, r in enumerate(fit_ranges)}
  for name, r in fit_ranges_dict.items():
    x.setRange(name, r[0], r[1])
  return ",".join(fit_ranges_dict.keys())

def robust_fit(extroopdf, pdf, datahist, fit_ranges_str, n_fits=8, recursive=True,
               max_n_fits=1024, seed=None):
  nlls = []
  free_params_vals = []
  
  log.info(f"Fitting {pdf.roopdf.GetName()} to {datahist.GetName()}. Doing {n_fits} fits from random initialisations.")
  for i in range(n_fits):
    pdf.randomize_params(None if seed is None else seed + i)
    r = pdf.roopdf.fitTo(datahist, Range=fit_ranges_str, PrintLevel=-1, Save=True, SumW2Error=True)
    nlls.append(r.minNll())
    free_params_vals.append(pdf.free_params_vals)

  best_free_params_vals = free_params_vals[np.argmin(nlls)]
  pdf.free_params_vals = best_free_params_vals

  max_diff = 0.01

  nll1 = min(nlls[:n_fits//2])
  nll2 = min(nlls[n_fits//2:])
  nll_diff = abs(nll1 - nll2)
  log.debug(f"Difference in minimum NLL from first and second half of fits is {nll_diff}")

  if nll_diff > max_diff:
    log.warning(f"Fit is unstable when starting from {n_fits} random initialisations.")
    
    if recursive and n_fits < max_n_fits:
      n_fits_more = n_fits * 2
      robust_fit(extroopdf, pdf, datahist, fit_ranges_str, n_fits_more, 
                 max_n_fits=max_n_fits, seed=seed)
    else:
      raise Exception(f"Fit is unstable. Tried {n_fits} from random places and did not find acceptable convergence. Halted because we reached maximum number of fits ({max_n_fits}).")

def fit(pdf, datahist, fit_ranges=(), method="robust", seed=None):
  fit_ranges_str = prepare_ranges(pdf.x, fit_ranges)

  # extended fit required to get valid results from fits in ranges (see https://root.cern/doc/v630/rf204b__extendedLikelihood__rangedFit_8py.html)
  n = ROOT.RooRealVar("n", "n", datahist.sumEntries(), 0, datahist.sumEntries())
  extroopdf = ROOT.RooAddPdf("extroopdf", "extroopdf", [pdf.roopdf], [n])

  assert method in ["robust", "from_defaults", "randomize"], f"Unknown fitting method: {method}"
  if method == "robust":
    robust_fit(extroopdf, pdf, datahist, fit_ranges_str, seed=seed)
  else:
    if method == "randomize":
      pdf.randomize_params(seed)
    r = extroopdf.fitTo(datahist, Range=fit_ranges_str, PrintLevel=-1, SumW2Error=True, Save=True)
    r.Print()

  pdf.check_bounds()
  
  twoNLL = 2*pdf.roopdf.createNLL(datahist, Range=fit_ranges_str, Offset="bin").getVal()
  
  fit_dof = int(utils.getNBinsFitted(pdf.x, fit_ranges) - pdf.get_dof())
  gof_pval = ROOT.TMath.Prob(twoNLL, fit_dof)
  return {"twoNLL": twoNLL, "gof_pval": gof_pval}

def main(in_file, out_file, pdf_name="Gaussian", order=1, fit_ranges=(), #
         nbins=None, plot_savepath=None, plot_range=None, method="robust"):
  log.info(f"Fitting {pdf_name} (order {order}) to events in {in_file} in ranges: {fit_ranges}")
  x, data = utils.readEvents(in_file)

  if nbins is None:
    nbins = x.getBins()
  if isinstance(data, ROOT.RooDataHist):
    assert x.getBins() % nbins == 0, f"Choosen number of bins ({nbins}) must be an exact divisor of the original binning (nbins = {x.getBins()})"
  x.setBins(nbins)

  log.info(f"Creating datahist with {nbins} bins")
  datahist = ROOT.RooDataHist("datahist", "datahist", ROOT.RooArgList(x), data)

  log.debug("Initialising fit function")
  pdf = getattr(pdfs, pdf_name)(x, postfix="cat0", order=order)
  fit(pdf, datahist, fit_ranges=fit_ranges, method=method)

  if plot_savepath is not None:
    xlim = (x.getMin(), x.getMax()) if plot_range == () else plot_range
    plotting.plotFit(datahist, pdf, plot_savepath, xlim=xlim)
    plotting.plotFitRoot(datahist, pdf, plot_savepath, xlim=xlim)

  log.info("Create output workspace")
  wsig = ROOT.RooWorkspace("wsig", "wsig")
  wsig.Import(datahist)
  wsig.Import(pdf.roopdf)
  log.info(f"Writing workspace to {out_file}")
  wsig.writeToFile(out_file)

if __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser(
                    prog='Signal Fitter',
                    description='Fits signal',
                    epilog='Text at the bottom of help')
  utils.addLoggingArguments(parser)
  parser.add_argument("in_file", type=str)
  parser.add_argument("out_file", type=str)
  parser.add_argument("--pdf-name", "-p", type=str, default="Gaussian", choices=pdfs.available_pdfs)
  parser.add_argument("--order", "-o", type=int, default=1, help="Function order")
  parser.add_argument("--plot-savepath", type=str, default=None)
  parser.add_argument("--plot-range", type=utils.comma_separated_two_tuple, default=(),
                      help="Range for plotting. Default is minimum and maximum value of x.")
  parser.add_argument("--fit-ranges", "-r", type=utils.comma_separated_two_tuple, nargs="+", default=())
  parser.add_argument("--nbins", type=int, default=None, help="Number of bins in to perform fit with. Default is the number of bins in input workspace.")
  parser.add_argument("--method", "-m", type=str, choices=["robust", "from_defaults", "randomize"], default="robust", help="")
  args = parser.parse_args()

  utils.applyLoggingArguments(args)  
  main(args.in_file, args.out_file, args.pdf_name, args.order, args.fit_ranges,
       args.nbins, args.plot_savepath, args.plot_range, args.method)