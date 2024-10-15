import os

import ROOT

def parse_path(path):
  ps = path.split(".root")
  assert len(ps) <= 2
  info = {
    "root_file_path": ps[0]+".root",
    "root_directory": os.path.dirname(ps[1]),
    "obj_path": os.path.basename(ps[1])
  }  
  return info

ws_obj_types = ["pdf", "var", "data", "function"]
obj_types = ["file", "workspace"] + ws_obj_types

def get_obj(info, obj_type):
  assert obj_type in obj_types
  
  f = ROOT.TFile(info["root_file_path"])
  f.cd(info["root_directory"])

  if obj_type == "file":
    return f

  ops = info["obj_path"].split(":")
  if len(ops) == 1 or obj_type == "workspace":
    return f.Get(ops[0])

  w = f.Get(ops[0])
  assert obj_type in ws_obj_types
  return getattr(w, obj_type)(ops[1])

def get_obj_from_path(path, obj_type):
  info = parse_path(path)
  return get_obj(info, obj_type)

def get_workspace(path):
  return get_obj_from_path(path, "workspace")

def get_file(path):
  return get_obj_from_path(path, "file")

def get_pdf(path):
  return get_obj_from_path(path, "pdf")

def get_pdf_with_x(path):
  info = parse_path(path)
  w = get_obj(info, "workspace")
  pdf_name = info["obj_path"].split(":")[1]
  x_name = "x"

  return w.pdf(pdf_name), w.var(x_name)

def get_var(path):
  return get_obj_from_path(path, "var")

def get_data(path):
  return get_obj_from_path(path, "data")

def get_function(path):
  return get_obj_from_path(path, "function")