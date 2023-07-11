import ROOT
import argparse
import common.functions as functions

def main(out_file, function="Gaussian", order=1, nevents=10000, randomise=False, n_datasets=0):
  x = ROOT.RooRealVar("x", "x", 100, 180)
  x.setBins(80)

  f = getattr(functions, function)(x, order=order, randomise=randomise)
  f.pdf.Print()
  for var in f.vars:
    var.Print()

  w = ROOT.RooWorkspace("w", "workspace")

  data = f.pdf.generateBinned(ROOT.RooArgSet(x), nevents)
  data.SetName("data")
  w.Import(data)

  for i in range(n_datasets):
    if i % 100 == 0:
      print(i)
    functions.randomiseVars(f.vars)
    data = f.pdf.generateBinned(ROOT.RooArgSet(x), nevents)
    data.SetName(f"data{i}")
    w.Import(data)

  w.writeToFile(out_file)

if __name__=="__main__":
  import argparse
  import common.args_and_logging as al
  parser = argparse.ArgumentParser(
                    prog='Signal Fit Plotter',
                    description='Plots signal fits')
  al.addLoggingArguments(parser)
  parser.add_argument("out_file")
  parser.add_argument("--function", "-f", type=str, default="Gaussian")
  parser.add_argument("--order", "-o", type=int, default=1)
  parser.add_argument("--nevents", "-n", type=int, default=10000)
  parser.add_argument("--randomise", action="store_true")
  parser.add_argument("--n-datasets", type=int, default=0)
  #parser.add_argument("--xlim", nargs=2, type=float, default=(115, 135))
  args = parser.parse_args()

  al.applyLoggingArguments(args)  
  main(args.out_file, args.function, args.order, args.nevents, args.randomise, args.n_datasets)