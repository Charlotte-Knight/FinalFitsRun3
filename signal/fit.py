import ROOT
import common.functions as functions
import common.tools as tools
import plotting.basic_fit
import logging
log = logging.getLogger(__name__)

def main(in_file, out_file, function="Gaussian", order=1, fit_ranges=[], plot_savepath=None):
  log.info(f"Will fit {function} (order {order}) to events in {in_file} in ranges: {fit_ranges}")
  x, dataset = tools.readEvents(in_file)

  for i, r in enumerate(fit_ranges):
    low, high = r.split(",")
    low, high = float(low), float(high)
    x.setRange(f"range{i}", low, high)
  x.setRange("Full", x.getMin(), x.getMax())
  fit_ranges_str = ",".join([f"range{i}" for i, r in enumerate(fit_ranges)])

  log.info("Creating datahist from dataset")
  datahist = ROOT.RooDataHist("datahist", "datahist", ROOT.RooArgList(x), dataset)

  log.info("Initialising fit function")
  f = getattr(functions, function)(x, postfix="cat0", order=order, randomise=True)
  log.info("Starting fit")
  N = ROOT.RooRealVar("N", "N", datahist.sumEntries(), 0, datahist.sumEntries()*2)
  extmodel = ROOT.RooExtendPdf("extmodel", "extmodel", f.pdf, N, "Full")
  extmodel.fitTo(datahist, ROOT.RooFit.Save(), ROOT.RooFit.Range(fit_ranges_str), ROOT.RooFit.PrintLevel(-1))

  if plot_savepath is not None:
    xlim = (x.getMin(), x.getMax())
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
  parser.add_argument("--function", "-f", type=str, default="Gaussian")
  parser.add_argument("--order", "-o", type=int, default=1)
  parser.add_argument("--plot-savepath","-p", type=str, default=None)
  parser.add_argument("--fit-ranges", type=str, nargs="+", default=[])
  args = parser.parse_args()

  print(args.in_file)

  al.applyLoggingArguments(args)  
  main(args.in_file, args.out_file, args.function, args.order, args.fit_ranges, args.plot_savepath)