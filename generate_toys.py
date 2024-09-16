import ROOT
import argparse
import common.functions as functions
from tqdm import tqdm
import logging
log = logging.getLogger(__name__)

def generateBinned(x, f, nevents, w, postfix="", randomise=False, asimov=False):
  if randomise:
    functions.randomiseVars(f.vars)
  data = f.pdf.generateBinned(x, nevents, ExpectedData=asimov)
  data.SetName(f"data{postfix}")
  w.Import(data)

def main(out_file, function="Gaussian", order=1, nevents=10000, randomise=False, 
         ndatasets=0, xlim=[100, 180], nbins=None, asimov=False):

  if nbins is None:
    nbins = int(xlim[1]-xlim[0]) # 1 bin per GeV

  log.info(f"Generating {'1' if ndatasets==0 else ndatasets}{' asimov' if asimov else ''} dataset(s) with {nevents} events from a {function} (order {order})")
  log.info(f"Dataset(s) have {nbins} bins between {xlim[0]} and {xlim[1]} GeV")
  log.info("Function parameters will be randomised" if randomise else "Function parameters will be kept to defaults")
  
  x = ROOT.RooRealVar("x", "x", xlim[0], xlim[1])
  x.setBins(nbins)

  f = getattr(functions, function)(x, order=order, randomise=randomise)
  log.debug(str(f.pdf).strip("\n"))
  for var in f.vars:
    log.debug(str(var).strip("\n"))

  w = ROOT.RooWorkspace("w", "workspace")

  if ndatasets == 0:
    generateBinned(x, f, nevents, w, randomise=randomise, asimov=asimov)
  else:
    for i in tqdm(range(ndatasets)):
      generateBinned(x, f, nevents, w, randomise=randomise, asimov=asimov, postfix=i)

  w.writeToFile(out_file)

if __name__=="__main__":
  import argparse
  import common.args_and_logging as al
  parser = argparse.ArgumentParser(prog='Binned Toy Generation')        
  al.addLoggingArguments(parser)
  parser.add_argument("out_file")
  parser.add_argument("--function", "-f", type=str, default="Gaussian", choices=functions.functions)
  parser.add_argument("--order", "-o", type=int, default=1, help="Function order")
  parser.add_argument("--nevents", "-n", type=int, default=10000, help="Number of events in a generated dataset")
  parser.add_argument("--randomise", action="store_true", help="Randomise the function parameters")
  parser.add_argument("--ndatasets", type=int, default=0, 
                      help="Set to non-zero value to generate multiple datasets. Will be named 'data0', 'data1'...")
  parser.add_argument("--xlim", nargs=2, type=float, default=(100,180), help="Limits on x. Default is 100,180.")
  parser.add_argument("--nbins", type=int, default=None, help="Number of bins in histogram. Default is 1/GeV.")
  parser.add_argument("--asimov", action="store_true", help="Generated asimov dataset(s)")
  args = parser.parse_args()

  al.applyLoggingArguments(args)  
  main(args.out_file, args.function, args.order, args.nevents, 
       args.randomise, args.ndatasets, args.xlim, args.nbins,
       args.asimov)