import logging

import numpy as np
import ROOT

from finalfits import pdfs, utils, fitting, plotting

log = logging.getLogger(__name__)

def get_roocategory(category_names):
  cat = ROOT.RooCategory("cat", "cat")
  for cat_name in category_names:
    cat.defineType(cat_name)
  return cat

def main(in_file, out_file, masses, pdf_name="Gaussian", order=1, plot_savepath=None):
  f = ROOT.TFile(in_file)
  w = f.Get("w")
  x = w.var("x")
  
  all_datahists = {f"{float(k.GetName().split('_')[-1]):.2f}": k for k in w.allData()}
  category_names = [f"{m:.2f}" for m in masses]
  assert np.isin(category_names, list(all_datahists)).all(), f"Masses in file: {list(all_datahists)}, required: {category_names}"
  
  datahists = {k: all_datahists[k] for k in category_names}
    
  x = ROOT.RooRealVar("x", "x", 110, 140)
  MH = ROOT.RooRealVar("MH", "MH", 125)
  MH.setConstant(True)
  pdf = getattr(pdfs, pdf_name)(x, order=order,
                                transforms={"mean*": [MH, 1]}, polys={"mean*|sigma*": [MH, 1]})
  
  cat = get_roocategory(category_names)
  w = ROOT.RooWorkspace("wtemp", "wtemp")
  w.Import({pdf.roopdf, cat})
  sct = ROOT.RooSimWSTool(w)
  simpdf = sct.build("simpdf", pdf.roopdf.GetName(), SplitParam=("MH", "cat"))
  
  for m in category_names:
    w.var(f"MH_{m}").setVal(float(m))
    
  combdata = ROOT.RooDataHist("combdata", "combdata", ROOT.RooArgList(x), Index=cat, Import=datahists)
  res = simpdf.fitTo(combdata, PrintLevel=-1, Save=True)
  res.Print()

  if plot_savepath:
    for m in category_names:
      w.var("MH").setVal(float(m))
      plotting.plotFit(datahists[m], w.pdf(pdf.roopdf.GetName()), w.var("x"), xlim=(110, 140))
    utils.savefig(plot_savepath)

  wout = ROOT.RooWorkspace("w", "w")
  wout.Import(w.pdf(pdf.roopdf.GetName()))
  wout.writeToFile(out_file)

if __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser(
                    prog='Signal Fitter',
                    description='Fits signal',
                    epilog='Text at the bottom of help')
  utils.addLoggingArguments(parser)
  parser.add_argument("in_file", type=str)
  parser.add_argument("out_file", type=str)
  parser.add_argument("--pdf-name", "-p", type=str, default="Gaussian", choices=["DCB", "Gaussian"])
  parser.add_argument("--order", "-o", type=int, default=1, help="Function order")
  parser.add_argument("--plot-savepath", type=str, default=None)
  parser.add_argument("--masses", "-m", nargs="+", type=float, default=[125], help="Masses to fit to")
  args = parser.parse_args()

  utils.applyLoggingArguments(args)  
  main(args.in_file, args.out_file, args.masses, args.pdf_name, args.order, args.plot_savepath)
