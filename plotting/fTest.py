import numpy as np
import matplotlib.pyplot as plt
import mplhep
mplhep.style.use("CMS")
import logging
log = logging.getLogger(__name__)

import common.tools as tools
import plotting.tools as ptools

def plotFamily(datahist, x, results, savepath, title=None):
  log.info("Plotting fTest fits")
  
  xlim = (100, 180)
  bin_width = ptools.histPlotTemplate(datahist, xlim)

  sf = datahist.sumEntries() * bin_width
  xi = np.linspace(xlim[0], xlim[1], 1000)
  for res in results:
    label = f"{title} {res['dof']}: " + r"$p_{ftest}=%.2f$, "%res["ftest_pval"] + r"$p_{gof}=%.2f$ "%res["gof_pval"]
    plt.plot(xi, tools.getVal(res["f"].pdf, x, xi)*sf, label=label)

  plt.legend()
  ptools.savefig(savepath)
  
def plotEnvelope(datahist, x, results, savepath):
  log.info("Plotting envelope")

  xlim = (100, 180)
  bin_width = ptools.histPlotTemplate(datahist, xlim)
  
  sf = datahist.sumEntries() * bin_width
  xi = np.linspace(xlim[0], xlim[1], 1000)

  gofs = [res["gof_pval"] for family in results.keys() for res in results[family]]
  best_gof_index = gofs.index(max(gofs))

  for family in results.keys():
    for res in results[family]:
      label = f"{family} {res['dof']}"
      plt.plot(xi, tools.getVal(res["f"].pdf, x, xi)*sf, label=label)

  legend = plt.legend()
  handle = legend.get_texts()[best_gof_index]
  handle.set_weight("bold")
  text = handle.get_text() + "\n(Best fit)"
  handle.set_text(text)

  ptools.savefig(savepath)