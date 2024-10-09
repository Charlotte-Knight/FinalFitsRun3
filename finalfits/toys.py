"""
Module to help with generating toy datasets. Primarily used for testing purposes.
"""

import logging

import ROOT
from tqdm import tqdm

from finalfits import pdfs, utils

log = logging.getLogger(__name__)

def generateBinned(x, pdf, nevents, w=None, postfix="", randomize=False, asimov=False):
  if randomize:
    pdf.randomize_params()
  data = pdf.roopdf.generateBinned(x, nevents, ExpectedData=asimov)
  data.SetName(f"data{postfix}")
  if w is not None:
    w.Import(data)
  return data

def main(out_file, pdf_name="Gaussian", order=1, nevents=10000, randomize=False,
         ndatasets=0, xlim=(100, 180), nbins=None, asimov=False):

  if nbins is None:
    nbins = int(xlim[1]-xlim[0]) # 1 bin per GeV

  log.info("Generating %s%s dataset(s) with %d events from a %s (order %d)",
           '1' if ndatasets==0 else ndatasets, ' asimov' if asimov else '',
           nevents, pdf_name, order)
  log.info("Dataset(s) have %d bins between %d and %d GeV", nbins, xlim[0], xlim[1])
  log.info("Function parameters will be %s", "randomized" if randomize else "kept to defaults")
    
  x = ROOT.RooRealVar("x", "x", xlim[0], xlim[1])
  x.setBins(nbins)

  pdf = getattr(pdfs, pdf_name)(x, order=order)
  if randomize:
    pdf.randomize_params()
  log.debug(str(pdf.roopdf).strip("\n"))
  for p in pdf.params:
    log.debug(str(p).strip("\n"))

  w = ROOT.RooWorkspace("w", "workspace")

  if ndatasets == 0:
    generateBinned(x, pdf, nevents, w, asimov=asimov)
  else:
    for i in tqdm(range(ndatasets)):
      generateBinned(x, pdf, nevents, w, randomize=randomize, asimov=asimov, postfix=i)

  w.writeToFile(out_file)

if __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser(prog='Binned Toy Generation')        
  utils.addLoggingArguments(parser)
  parser.add_argument("out_file")
  parser.add_argument("--pdf-name", "-p", type=str, default="Gaussian", choices=pdfs.available_pdfs)
  parser.add_argument("--order", "-o", type=int, default=1, help="Function order")
  parser.add_argument("--nevents", "-n", type=int, default=10000, help="Number of events in a generated dataset")
  parser.add_argument("--randomize", action="store_true", help="randomize the function parameters")
  parser.add_argument("--ndatasets", type=int, default=0, 
                      help="Set to non-zero value to generate multiple datasets. Will be named 'data0', 'data1'...")
  parser.add_argument("--xlim", type=utils.comma_separated_two_tuple, default=(100,180), help="Limits on x. Default is 100,180.")
  parser.add_argument("--nbins", type=int, default=None, help="Number of bins in histogram. Default is 1/GeV.")
  parser.add_argument("--asimov", action="store_true", help="Generated asimov dataset(s)")
  args = parser.parse_args()

  utils.applyLoggingArguments(args)  
  main(args.out_file, args.pdf_name, args.order, args.nevents,
       args.randomize, args.ndatasets, args.xlim, args.nbins,
       args.asimov)