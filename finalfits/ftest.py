import logging

import ROOT

from finalfits import plotting, fitting, utils, pdfs

log = logging.getLogger(__name__)

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

def getResults(x, datahist, pdf_name, max_dof, fit_ranges, blinded_regions, plot_savepath, do_all_orders):
  results = []
  order = 1
  while shouldKeepGoing(results, do_all_orders, max_dof):
    pdf = getattr(pdfs, pdf_name)(x, postfix="cat0", order=order)
    
    pdf_info = {"pdf": pdf, "order": pdf.order, "dof": pdf.get_dof()}
    fit_result = fitting.fit(pdf, datahist, fit_ranges)
    results.append(pdf_info | fit_result)    
    
    log.info(f"Goodness-of-fit pval = {results[-1]['gof_pval']:.2f}")

    if len(results) == 1:
      ftest_pval = 0
    else:
      dchi2 = results[-2]["twoNLL"]-results[-1]["twoNLL"]
      if dchi2 < 0:
        dchi2 = 0
      log.debug(f"Delta Chi2 = {dchi2:.2f}")
      ftest_pval = ROOT.TMath.Prob(dchi2, results[-1]["dof"]-results[-2]["dof"])
      
    results[-1]["ftest_pval"] = ftest_pval
    log.info(f"F-test pval = {ftest_pval:.2f}")

    order += 1
  
  if plot_savepath is not None:
    plotting.plotFamily(datahist, x, results, plot_savepath, blinded_regions, pdf_name)

  return results

def filterByGof(results, gof_threshold=0.01):
  new_results = {}
  for family in results.keys():
    new_results[family] = [res for res in results[family] if res["gof_pval"] > gof_threshold]
    if len(new_results[family]) > 0:
      new_results[family][0]["ftest_pval"] = 0.0 # must accept first pdf that passes gof
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
  pdfs = ROOT.RooArgList(*[res["pdf"].roopdf for res in flattened_results])
  gofs = [res["gof_pval"] for family in results.keys() for res in results[family]]
  
  pdfIndex = ROOT.RooCategory("pdfIndex", "pdfIndex")
  multipdf = ROOT.RooMultiPdf("multipdf", "multipdf", pdfIndex, pdfs)
  pdfIndex.setIndex(gofs.index(max(gofs)))

  return pdfIndex, multipdf

def main(in_file, out_file, pdf_names, max_dof=5, fit_ranges=[], blinded_regions=[], plot_savepath=None, do_all_orders=False):
  x, dataset = utils.readEvents(in_file)
  
  log.info("Creating datahist from dataset")
  datahist = ROOT.RooDataHist("datahist", "datahist", ROOT.RooArgList(x), dataset)

  results = {pdf_name: getResults(x, datahist, pdf_name, max_dof, fit_ranges, blinded_regions, None if plot_savepath is None else plot_savepath+pdf_name, do_all_orders) for pdf_name in pdf_names}
  results = filterResults(results)
  if plot_savepath is not None:
    plotting.plotEnvelope(datahist, x, results, plot_savepath+"Envelope", blinded_regions)

  createEnvelope(results)

if __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser(
                    prog='Signal Fitter',
                    description='Fits signal',
                    epilog='Text at the bottom of help')
  utils.addLoggingArguments(parser)
  parser.add_argument("in_file", type=str)
  parser.add_argument("out_file", type=str)
  parser.add_argument("--pdf-names", "-p", type=str, nargs="+", default=["Power", "Exponential", "ExpPoly", "Bernstein", "Laurent"], choices=pdfs.available_pdfs)
  parser.add_argument("--max-dof", "-o", type=int, default=5)
  parser.add_argument("--plot-savepath", type=str, default=None)
  parser.add_argument("--fit-ranges", type=utils.comma_separated_two_tuple, nargs="+", default=[(100,120), (130,180)])
  parser.add_argument("--blinded-regions", type=str, nargs="+", default=["115,135"])
  parser.add_argument("--do-all-orders", action="store_true")
  args = parser.parse_args()

  utils.applyLoggingArguments(args)  
  main(args.in_file, args.out_file, args.pdf_names, args.max_dof, args.fit_ranges, 
       args.blinded_regions, args.plot_savepath, args.do_all_orders)