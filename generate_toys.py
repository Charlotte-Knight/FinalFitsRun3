import ROOT
import argparse
import common.functions as functions
from tqdm import tqdm
import common.tools as tools

def generateBinned(x, f, nevents, w, postfix="", randomise=False):
  if randomise:
    functions.randomiseVars(f.vars)
  #data = f.pdf.generateBinned(x, nevents)
  data = tools.makeAsimovDataHist(x, f.pdf, norm=nevents)
  data.SetName(f"data{postfix}")
  w.Import(data)

def main(out_file, function="Gaussian", order=1, nevents=10000, randomise=False, 
         ndatasets=0, xlim=[100, 180], nbins=None):

  if nbins is None:
    nbins = int(xlim[1]-xlim[0])

  x = ROOT.RooRealVar("x", "x", xlim[0], xlim[1])
  x.setBins(nbins)

  f = getattr(functions, function)(x, order=order, randomise=randomise)
  f.pdf.Print()
  for var in f.vars:
    var.Print()

  w = ROOT.RooWorkspace("w", "workspace")

  generateBinned(x, f, nevents, w, randomise=randomise)
  for i in tqdm(range(ndatasets)):
    generateBinned(x, f, nevents, w, randomise=randomise, postfix=i)

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
  parser.add_argument("--ndatasets", type=int, default=0)
  parser.add_argument("--xlim", nargs=2, type=float, default=(100, 180))
  parser.add_argument("--nbins", type=int, default=None)
  args = parser.parse_args()

  al.applyLoggingArguments(args)  
  main(args.out_file, args.function, args.order, args.nevents, 
       args.randomise, args.ndatasets, args.xlim, args.nbins)