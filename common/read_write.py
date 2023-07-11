import ROOT

def parsePath(path):
  # default values
  info = {"root_directory": "", "workspace_name":"w", "data_name":"data"}

  root_file_split = path.split(".root")
  assert len(root_file_split) <= 2
  info["root_file_path"] = root_file_split[0]+".root"

  if len(root_file_split) == 1:
    return info

  info["root_directory"] = "/".join(root_file_split[1].split("/")[:-1]) + "/"
  info["workspace_name"] = root_file_split[1].split("/")[-1].split(":")[0]
  info["data_name"] = root_file_split[1].split(":")[1]
  
  return info

def getRootFile(path):
  info = parsePath(path)
  f = ROOT.TFile(info["root_file_path"])
  return f

def getWorkspace(path):
  info = parsePath(path)
  f = ROOT.TFile(info["root_file_path"])
  assert f.cd(info["root_directory"])
  w = f.Get(info["workspace_name"])
  return w

def getData(path):
  info = parsePath(path)
  f = ROOT.TFile(info["root_file_path"])
  assert f.cd(info["root_directory"])
  w = f.Get(info["workspace_name"])
  data = w.data("")
  return w