import ROOT
import logging

finalfits_verbose_dict = {
  -2: "CRITICAL",
  -1: "ERROR",
  0:  "WARNING",
  1:  "INFO",
  2:  "DEBUG"
}

roofit_verbose_dict = {
  -2: "FATAL",
  -1: "ERROR",
  0:  "WARNING",
  1:  "PROGRESS",
  2:  "INFO",
  3:  "DEBUG"
}

def addLoggingArguments(parser):
  parser.add_argument("--verbose", "-v", type=int, default=1, choices=range(-2,3),
                      help="Set verbosity level for finalFits scripts: %s"%(", ".join(f"{key}={value}" for key,value in finalfits_verbose_dict.items())))
  parser.add_argument("--roofit-verbose", "-rv", type=int, default=0, choices=range(-2,4),
                      help="Set verbosity level for RooFit: %s"%(", ".join(f"{key}={value}" for key,value in roofit_verbose_dict.items())))

def applyLoggingArguments(args):
  ROOT.RooMsgService.instance().setGlobalKillBelow(getattr(ROOT.RooFit, roofit_verbose_dict[args.roofit_verbose]))
  logging.basicConfig(level=getattr(logging, finalfits_verbose_dict[args.verbose]), format=('%(name)-20s: %(levelname)-8s %(message)s'))
