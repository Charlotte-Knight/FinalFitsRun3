import logging

import ROOT

from finalfits import pdfs, utils, toys

log = logging.getLogger(__name__)

def main(out_file, masses, pdf_name="Gaussian", order=1, nevents=10000, randomize=False,
         xlim=(100, 180), nbins=None, asimov=False):

  if nbins is None:
    nbins = int(xlim[1]-xlim[0]) # 1 bin per GeV

  log.info("Generating%s dataset(s) for masses: %s with %d events from a %s (order %d)",
           ' asimov' if asimov else '', masses, nevents, pdf_name, order)
  log.info("Dataset(s) have %d bins between %d and %d GeV", nbins, xlim[0], xlim[1])
  log.info("Function parameters will be %s", "randomized" if randomize else "kept to defaults")
    
  x = ROOT.RooRealVar("x", "x", xlim[0], xlim[1])
  x.setBins(nbins)

  MH = ROOT.RooRealVar("MH", "MH", 125)
  MH.setConstant(True)
  pdf = getattr(pdfs, pdf_name)(x, order=order, transforms={"mean.*": [MH, 1]}, 
                                polys={"mean.*": [MH, 1], "sigma.*": [MH, 1]})
  if randomize:
    pdf.randomize_params()
  log.debug(str(pdf.roopdf).strip("\n"))
  for p in pdf.params:
    log.debug(str(p).strip("\n"))

  w = ROOT.RooWorkspace("w", "workspace")

  for m in masses:
    pdf["MH"] = m
    toys.generateBinned(x, pdf, nevents, w, asimov=asimov, postfix=f"_{m}")

  w.writeToFile(out_file)

if __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser(prog='Binned Toy Generation')        
  utils.addLoggingArguments(parser)
  parser.add_argument("out_file")
  parser.add_argument("--pdf-name", "-p", type=str, default="Gaussian", choices=["DCB", "Gaussian"])
  parser.add_argument("--order", "-o", type=int, default=1, help="Function order")
  parser.add_argument("--nevents", "-n", type=int, default=10000, help="Number of events in a generated dataset")
  parser.add_argument("--randomize", action="store_true", help="randomize the function parameters")
  parser.add_argument("--xlim", type=utils.comma_separated_two_tuple, default=(100,180), help="Limits on x. Default is 100,180.")
  parser.add_argument("--nbins", type=int, default=None, help="Number of bins in histogram. Default is 1/GeV.")
  parser.add_argument("--asimov", action="store_true", help="Generated asimov dataset(s)")
  parser.add_argument("--masses", "-m", nargs="+", type=float, default=[125], help="Masses to generate datasets for")
  args = parser.parse_args()
 
  utils.applyLoggingArguments(args)  
  main(args.out_file, args.masses, args.pdf_name, args.order, args.nevents,
       args.randomize, args.xlim, args.nbins, args.asimov)