import argparse

from finalfits import utils, plotting

parser = argparse.ArgumentParser(prog='Binned Toy Generation')        
utils.addLoggingArguments(parser)
parser.add_argument("in_file")
parser.add_argument("out_file")
args = parser.parse_args()
utils.applyLoggingArguments(args)  

x, data = utils.readEvents(args.in_file)
plotting.plotHist(data, x, args.out_file)
