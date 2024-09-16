import ROOT
import common.functions as functions
import common.tools as tools
import plotting.basic_fit
import logging
log = logging.getLogger(__name__)

def main(in_file, out_file, function="Gaussian", order=1, fit_ranges=[], #
         nbins=None, plot_savepath=None, plot_range=None):
  log.info(f"Fitting {function} (order {order}) to events in {in_file} in ranges: {fit_ranges}")
  x, data = tools.readEvents(in_file)

  if nbins is None:
    nbins = x.getBins()
  if isinstance(data, ROOT.RooDataHist):
    assert x.getBins() % nbins == 0, f"Choosen number of bins ({nbins}) must be an exact divisor of the original binning (nbins = {x.getBins()})"
  x.setBins(nbins)

  # prepare fit ranges for RooFit
  for i, r in enumerate(fit_ranges):
    low, high = r.split(",")
    low, high = float(low), float(high)
    x.setRange(f"range{i}", low, high)
  x.setRange("Full", x.getMin(), x.getMax())
  fit_ranges_str = ",".join([f"range{i}" for i, r in enumerate(fit_ranges)]) # to be passed to fitTo

  log.info(f"Creating datahist with {nbins} bins")
  datahist = ROOT.RooDataHist("datahist", "datahist", ROOT.RooArgList(x), data)

  log.debug("Initialising fit function")
  f = getattr(functions, function)(x, postfix="cat0", order=order, randomise=True)
  N = ROOT.RooRealVar("N", "N", datahist.sumEntries(), 0, datahist.sumEntries()*2)
  extmodel = ROOT.RooExtendPdf("extmodel", "extmodel", f.pdf, N, "Full")

  log.info("Starting fit")
  extmodel.fitTo(datahist, SumW2Error=True, Range=fit_ranges_str, Save=True, PrintLevel=-1)

  if plot_savepath is not None:
    xlim = (x.getMin(), x.getMax()) if args.plot_range is None else args.plot_range
    plotting.basic_fit.plotFit(datahist, x, f.pdf, f.vars, plot_savepath, xlim=xlim)
    plotting.basic_fit.plotFitRoot(datahist, x, f.pdf, f.vars, plot_savepath, xlim=xlim)

  log.info("Create output workspace")
  wsig = ROOT.RooWorkspace("wsig", "wsig")
  wsig.Import(datahist)
  wsig.Import(f.pdf)
  log.info(f"Writing workspace to {out_file}")
  wsig.writeToFile(out_file)

if __name__=="__main__":
  import argparse
  import common.args_and_logging as al
  parser = argparse.ArgumentParser(
                    prog='Signal Fitter',
                    description='Fits signal',
                    epilog='Text at the bottom of help')
  al.addLoggingArguments(parser)
  parser.add_argument("in_file", type=str)
  parser.add_argument("out_file", type=str)
  parser.add_argument("--function", "-f", type=str, default="Gaussian", choices=functions.functions)
  parser.add_argument("--order", "-o", type=int, default=1, help="Function order")
  parser.add_argument("--plot-savepath","-p", type=str, default=None)
  parser.add_argument("--plot-range", nargs=2, type=float, default=None, 
                      help="Range for plotting. Default is minimum and maximum value of x.")
  parser.add_argument("--fit-ranges", type=str, nargs="+", default=[])
  parser.add_argument("--nbins", type=int, default=None, help="Number of bins in to perform fit with. Default is the number of bins in input workspace.")
  args = parser.parse_args()

  al.applyLoggingArguments(args)  
  main(args.in_file, args.out_file, args.function, args.order, args.fit_ranges, 
       args.nbins, args.plot_savepath, args.plot_range)