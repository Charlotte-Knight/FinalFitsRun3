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

def plot_function(pdf, function, fit_results, plot_savepath):
  single_fit_results = {k: fit_results[k] for k in fit_results if k != "sim"}
  MHs = [float(k) for k in single_fit_results]
  
  values = []
  errors = []
  for m, r in single_fit_results.items():
    pdf.w.loadSnapshot(f"fit_{m}")
    values.append(pdf[function].getVal())
    errors.append(pdf[function].getPropagatedError(r))
  hist = ROOT.TH1F("hist", "Histogram with errors", len(MHs), MHs[0], MHs[-1])
  
  print(values)
  
  for i, MH in enumerate(MHs):
    hist.SetBinContent(i, values[i])
    hist.SetBinError(i, errors[i])  
  data_hist = ROOT.RooDataHist("data_hist", "RooDataHist from TH1", [pdf["MH"]], hist)
  
  frame = pdf["MH"].frame()
  
  pdf.w.loadSnapshot("fit_sim")
  #pdf[function].plotOn(frame, VisualizeError=(fit_results["sim"], 1), FillColor="kOrange")
  #pdf[function].plotOn(frame)
  
  data_hist.plotOn(frame)


  c = ROOT.TCanvas("rf610_visualerror", "rf610_visualerror", 800, 800)
  frame.Draw()
  c.SaveAs(f"{plot_savepath}_{function}.png")

def main(in_file, out_file, masses, pdf_name="Gaussian", order=1, plot_savepath=None):
  f = ROOT.TFile(in_file)
  w = f.Get("w")
  x = w.var("x")
  
  all_datahists = {f"{float(k.GetName().split('_')[-1]):.2f}": k for k in w.allData()}
  category_names = [f"{m:.2f}" for m in masses]
  assert np.isin(category_names, list(all_datahists)).all(), f"Masses in file: {list(all_datahists)}, required: {category_names}"
  
  datahists = {k: all_datahists[k] for k in category_names}
    
  x = ROOT.RooRealVar("x", "x", 110, 140)
  MH = ROOT.RooRealVar("MH", "MH", 125, 100, 180)
  MH.setConstant(True)
  pdf = getattr(pdfs, pdf_name)(x, order=order, transforms={"mean.*": [MH, 1]},
                                polys={"mean.*": [MH, 1], "sigma.*": [MH, 1]})
  
  # single mass fits
  pdf.freeze_polys()
  fit_results = {}
  for m in category_names:
    pdf["MH"] = float(m)
    fit_results[m] = pdf.extroopdf.fitTo(datahists[m], PrintLevel=-1, Save=True, SumW2Error=True, Offset=True)
    fit_results[m].Print()
    pdf.w.saveSnapshot(f"fit_{m}", pdf.w.allVars())
  pdf.freeze_polys(False)
  
  pdf.w.Print()
  
  # build simultaneous pdf
  cat = get_roocategory(category_names)
  pdf.w.Import(cat)
  sct = ROOT.RooSimWSTool(pdf.w)
  simpdf = sct.build("simpdf", pdf.extroopdf.GetName(), SplitParam=("MH", "cat"))
  for m in category_names:
    pdf[f"MH_{m}"] = float(m)
  
  # fit simultaneous pdf
  combdata = ROOT.RooDataHist("combdata", "combdata", [x], Index=cat, Import=datahists)
  fit_results["sim"] = pdf.fit_result = simpdf.fitTo(combdata, PrintLevel=-1, Save=True, SumW2Error=True, Offset=True)
  pdf.w.saveSnapshot("fit_sim", pdf.w.allVars())
  
  if plot_savepath:
    for m in category_names:
      pdf["MH"] = float(m)
      plotting.plotFit(datahists[m], pdf, pdf.x, xlim=(110, 140))
    utils.savefig(plot_savepath)

  # plot functions  
  #plot_function(pdf, "mean1", fit_results, plot_savepath)
  #plot_function(pdf, "mean1_free", fit_results, plot_savepath)
  plot_function(pdf, "sigmaLR1", fit_results, plot_savepath)

  # wout = ROOT.RooWorkspace("w", "w")
  # wout.Import(w.pdf(pdf.roopdf.GetName()))
  # wout.writeToFile(out_file)

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
