import ROOT
import common.functions as functions
import plotting.signal_fit

def readEvents(filename):
  f = ROOT.TFile(filename)
  w = f.Get("w")
  x = w.var("x")
  dataset = w.data("data")
  return x, dataset

def main(in_file, out_file, function="Gaussian", order=1, plot_savepath=None):
  x, dataset = readEvents(in_file)
  datahist = ROOT.RooDataHist("datahist", "datahist", ROOT.RooArgList(x), dataset)

  f = getattr(functions, function)
  pdf, others = f(x, postfix="cat0", order=order)
  params = pdf.getParameters(x)
  pdf.fitTo(datahist, ROOT.RooFit.Save())

  if plot_savepath is not None:
    plotting.signal_fit.plotFit(datahist, x, pdf, params, plot_savepath, xlim=(110, 140))

  wsig = ROOT.RooWorkspace("wsig", "wsig")
  wsig.Import(datahist)
  wsig.Import(pdf)
  wsig.writeToFile(out_file)

if __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser(
                    prog='Signal Fitter',
                    description='Fits signal',
                    epilog='Text at the bottom of help')
  parser.add_argument("in_file")
  parser.add_argument("out_file")
  parser.add_argument("--function", "-f", type=str, default="Gaussian")
  parser.add_argument("--order", type=int, default=1)
  parser.add_argument("--plot-savepath","-p", type=str, default=None)
  args = parser.parse_args()

  main(args.in_file, args.out_file, args.function, args.order, args.plot_savepath)