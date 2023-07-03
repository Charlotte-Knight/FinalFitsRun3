import ROOT
import common.functions as functions
import plotting.signal_fit

def readEvents(filename):
  f = ROOT.TFile(filename)
  w = f.Get("w")
  x = w.var("x")
  dataset = w.data("data")
  return x, dataset

def main(in_file, out_file, function="Gaussian", order=1, fit_ranges=[], plot_savepath=None):
  x, dataset = readEvents(in_file)

  x.setRange("Full", x.getMin(), x.getMax())
  for i, r in enumerate(fit_ranges):
    low, high = r.split(",")
    low, high = float(low), float(high)
    x.setRange(f"range{i}", low, high)
  fit_ranges_str = ",".join([f"range{i}" for i, r in enumerate(fit_ranges)])

  datahist = ROOT.RooDataHist("datahist", "datahist", ROOT.RooArgList(x), dataset)

  f = getattr(functions, function)(x, postfix="cat0", order=order, randomise=True)
  f.pdf.Print()
  f.pdf.selectNormalizationRange("Full", True)
  #f.pdf.selectNormalizationRange(fit_ranges_str, True)
  res = f.pdf.fitTo(datahist, ROOT.RooFit.Save(), ROOT.RooFit.Range(fit_ranges_str))
  res.Print()

  if plot_savepath is not None:
    xlim = (x.getMin(), x.getMax())
    #xlim = (115, 135)
    plotting.signal_fit.plotFit(datahist, x, f.pdf, f.vars, plot_savepath, xlim=xlim)

  wsig = ROOT.RooWorkspace("wsig", "wsig")
  wsig.Import(datahist)
  wsig.Import(f.pdf)
  wsig.writeToFile(out_file)
  wsig.Print()

if __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser(
                    prog='Signal Fitter',
                    description='Fits signal',
                    epilog='Text at the bottom of help')
  parser.add_argument("in_file")
  parser.add_argument("out_file")
  parser.add_argument("--function", "-f", type=str, default="Gaussian")
  parser.add_argument("--order", "-o", type=int, default=1)
  parser.add_argument("--plot-savepath","-p", type=str, default=None)
  parser.add_argument("--fit-ranges", type=str, nargs="+", default=[])
  args = parser.parse_args()

  main(args.in_file, args.out_file, args.function, args.order, args.fit_ranges, args.plot_savepath)