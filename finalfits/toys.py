"""
Module to help with generating toy datasets. Primarily used for testing purposes.
"""

import logging

log = logging.getLogger(__name__)

def generateBinned(x, pdf, nevents, w=None, postfix="", randomize=False, asimov=False):
  if randomize:
    pdf.randomize_params()
  data = pdf.roopdf.generateBinned(x, nevents, ExpectedData=asimov)
  data.SetName(f"data{postfix}")
  if w is not None:
    w.Import(data)
  return data