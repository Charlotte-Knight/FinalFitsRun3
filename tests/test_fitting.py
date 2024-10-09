import pytest

import ROOT

from finalfits import pdfs, toys, fitting, plotting

fitting_tests = [
  (pdf_name, order, method, "unblinded", 0.02)
  for order in range(1, 6) 
  for pdf_name in pdfs.available_pdfs 
  for method in ["robust", "from_defaults", "randomize"]
  if order <= getattr(pdfs, pdf_name).max_order
]

fitting_tests += [
  (pdf_name, order, method, "blinded", 0.05)
  for order in range(1, 2)
  for pdf_name in pdfs.available_pdfs
  for method in ["robust", "from_defaults", "randomize"]
  if order <= getattr(pdfs, pdf_name).max_order
]

@pytest.mark.parametrize("pdf_name,order,method,blind_status,chi2_threshold", fitting_tests)
def test_fit(pdf_name, order, method, blind_status, chi2_threshold):
  r = (110, 140) if pdf_name in ["Gaussian", "DCB"] else (100, 180)
  nbins = 80
  x = ROOT.RooRealVar("x", "x", r[0], r[1])
  x.setBins(nbins)

  if blind_status == "blinded":
    # blind in the middle by 20% of the range
    mid = (r[0] + r[1]) / 2
    size = (r[1] - r[0]) / 5
    fit_ranges = ((r[0], mid-size), (mid+size, r[1]))
  else:
    fit_ranges = ((r[0], r[1]), )

  pdf = getattr(pdfs, pdf_name)(x, order=order)
  pdf.randomize_params(seed=0)
  datahist = toys.generateBinned(x, pdf, 100000, asimov=True)

  fitting.fit(pdf, datahist, method=method, fit_ranges=fit_ranges)
  
  chi2 = 2*pdf.roopdf.createChi2(datahist).getVal()
  dof = int(nbins - pdf.get_dof())
  chi2_dof = chi2 / dof
  
  if chi2_dof > chi2_threshold:
    plotting.plotFit(datahist, x, pdf.roopdf, pdf.params, f"tests/plots/{pdf_name}{order}_{method}_{blind_status}")
    
  assert chi2_dof <= chi2_threshold