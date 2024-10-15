import logging
import argparse

from finalfits import utils, plotting, read_write as rw

log = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
                  prog='Signal Fitter',
                  description='Fits signal',
                  epilog='Text at the bottom of help')
utils.addLoggingArguments(parser)
parser.add_argument("data_path", type=str)
parser.add_argument("pdf_path", type=str)
parser.add_argument("plot_savepath", type=str)
args = parser.parse_args()

utils.applyLoggingArguments(args)  

datahist = rw.get_data(args.data_path)
w = rw.get_workspace(args.pdf_path)
pdf = w.pdf("DCB1")
x = w.var("x")
MH = w.var("MH")
MH.setVal(126)

plotting.plotFit(datahist, pdf, x, args.plot_savepath)
