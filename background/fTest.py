import ROOT
import common.functions as functions
import plotting.fTest
import numpy as np
import logging
import common.tools as tools
log = logging.getLogger(__name__)

def fitFunction(x, datahist, function, order, fit_ranges_str, n_bins_fitted):
  log.info(f"Initialising {function} (order {order})")
  f = getattr(functions, function)(x, postfix="cat0", order=order)

  tools.robustFit(f.pdf, datahist, fit_ranges_str)
  
  twoNLL = 2*f.pdf.createNLL(datahist, Range=fit_ranges_str, Offset="bin").getVal()
  dof = int(n_bins_fitted - f.getDof())
  gof = ROOT.TMath.Prob(twoNLL, dof)

  return {"f":f, "order": f.order, "dof": f.getDof(), "twoNLL": twoNLL, "gof_pval": gof, "ftest_pval":0}

def shouldStop(results, gof_threshold=0.01, ftest_threshold=0.05):
  return (results[-1]["ftest_pval"] > ftest_threshold) and (results[-1]["gof_pval"] > gof_threshold)

def shouldKeepGoing(results, do_all_orders, max_dof, gof_threshold=0.01, ftest_threshold=0.05):
  if results == []:
    return True
  elif results[-1]["dof"] == max_dof:
    return False
  elif results[-1]["dof"] > max_dof:
    del results[-1]
    return False
  elif (not do_all_orders) and (results[-1]["ftest_pval"] > ftest_threshold) and (results[-1]["gof_pval"] > gof_threshold):
    return False
  else:
    return True

def getResults(x, datahist, function, max_dof, fit_ranges_str, blinded_regions, n_bins_fitted, plot_savepath, do_all_orders):
  results = []
  order = 1
  while shouldKeepGoing(results, do_all_orders, max_dof):
    results.append(fitFunction(x, datahist, function, order, fit_ranges_str, n_bins_fitted))
    log.info(f"Goodness-of-fit pval = {results[-1]['gof_pval']:.2f}")

    if len(results) > 1:
      dchi2 = results[-2]["twoNLL"]-results[-1]["twoNLL"]
      if dchi2 < 0:
        dchi2 = 0
      log.debug(f"Delta Chi2 = {dchi2:.2f}")
      ftest_pval = ROOT.TMath.Prob(dchi2, results[-1]["order"]-results[-2]["order"])
      results[-1]["ftest_pval"] = ftest_pval
      log.info(f"F-test pval = {ftest_pval:.2f}")

    order += 1
  
  if plot_savepath is not None:
    plotting.fTest.plotFamily(datahist, x, results, plot_savepath, blinded_regions, function)

  return results

def filterByGof(results, gof_threshold=0.01):
  new_results = {}
  for family in results.keys():
    new_results[family] = [res for res in results[family] if res["gof_pval"] > gof_threshold]
    if len(new_results[family]) > 0:
      new_results[family][0]["ftest_pval"] = 0.0 # must accept first function that passes gof
  return new_results

def filterByFtest(results, ftest_threshold=0.05):
  new_results = {}
  for family in results.keys():
    new_results[family] = [res for res in results[family] if res["ftest_pval"] < ftest_threshold]
  return new_results

def filterResults(results, gof_threshold=0.01, ftest_threshold=0.05):
  results = filterByGof(results, gof_threshold)
  results = filterByFtest(results, ftest_threshold)
  return results

def createEnvelope(results):
  flattened_results = [res for family in results.keys() for res in results[family]]
  pdfs = ROOT.RooArgList(*[res["f"].pdf for res in flattened_results])
  gofs = [res["gof_pval"] for family in results.keys() for res in results[family]]
  
  pdfIndex = ROOT.RooCategory("pdfIndex", "pdfIndex")
  multipdf = ROOT.RooMultiPdf("multipdf", "multipdf", pdfIndex, pdfs)
  pdfIndex.setIndex(gofs.index(max(gofs)))

  return pdfIndex, multipdf

def main(in_file, out_file, functions, max_dof=5, fit_ranges=[], blinded_regions=[], plot_savepath=None, do_all_orders=False):
  x, dataset = tools.readEvents(in_file)
  fit_ranges_str = tools.setRanges(x, fit_ranges)
  n_bins_fitted = tools.getNBinsFitted(x, fit_ranges)

  log.info("Creating datahist from dataset")
  datahist = ROOT.RooDataHist("datahist", "datahist", ROOT.RooArgList(x), dataset)

  results = {function: getResults(x, datahist, function, max_dof, fit_ranges_str, blinded_regions, n_bins_fitted, None if plot_savepath is None else plot_savepath+function, do_all_orders) for function in functions}
  results = filterResults(results)
  if plot_savepath is not None:
    plotting.fTest.plotEnvelope(datahist, x, results, plot_savepath+"Envelope", blinded_regions)

  createEnvelope(results)

if __name__=="__main__":
  import argparse
  import common.args_and_logging as al
  parser = argparse.ArgumentParser(
                    prog='Signal Fitter',
                    description='Fits signal',
                    epilog='Text at the bottom of help')
  al.addLoggingArguments(parser)
  parser.add_argument("in_file", type=str)
  parser.add_argument("out_file", type=str)
  parser.add_argument("--functions", "-f", type=str, nargs="+", default=["Power", "Exponential", "ExpPoly", "Bernstein"])
  parser.add_argument("--max-dof", "-o", type=int, default=5)
  parser.add_argument("--plot-savepath","-p", type=str, default=None)
  parser.add_argument("--fit-ranges", type=str, nargs="+", default=["100,115", "135,180"])
  parser.add_argument("--blinded-regions", type=str, nargs="+", default=["115,135"])
  parser.add_argument("--do-all-orders", action="store_true")
  args = parser.parse_args()

  al.applyLoggingArguments(args)  
  main(args.in_file, args.out_file, args.functions, args.max_dof, args.fit_ranges, 
       args.blinded_regions, args.plot_savepath, args.do_all_orders)