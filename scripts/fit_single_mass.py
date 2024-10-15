import logging

import ROOT
from finalfits import pdfs, utils, fitting, plotting

log = logging.getLogger(__name__)

def main(in_file, out_file, pdf_name="Gaussian", order=1, plot_savepath=None):
  x, datahist = utils.readEvents(in_file)

  MH = ROOT.RooRealVar("MH", "MH", 125)
  MH.setConstant(True)

  transforms = {"mean*": [MH, 1]}

  pdf = getattr(pdfs, pdf_name)(x, order=order, transforms=transforms)
  res = fitting.fit(pdf, datahist, fit_ranges=[(115, 135)], method="robust")
  log.info(f"Fit results: {res}")

  if plot_savepath:
    plotting.plotFit(datahist, pdf, pdf.x, plot_savepath, xlim=(115, 135))

  wsig = ROOT.RooWorkspace("w", "w")
  wsig.Import(datahist)
  wsig.Import(pdf.roopdf)
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
  args = parser.parse_args()

  utils.applyLoggingArguments(args)  
  main(args.in_file, args.out_file, args.pdf_name, args.order, args.plot_savepath)
