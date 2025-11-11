import matplotlib.pyplot as plt
import numpy as np

def load_file(fname: str):
  
  with open("/gluster/data/atlas/jdombrowski/SelectionCutIntegrals/" + fname + ".txt", "r") as file:
    file_string = file.read().split("\n");
  n = len(file_string)-1
  data = [[0 for i in range(1, n)], [0 for i in range(1, n)]]
  for i in range(1, n):
    data[0][i-1] = float(file_string[i].split(",")[0])
    data[1][i-1] = float(file_string[i].split(",")[1])
  return data

def load_dict(data_dict: dict):
  for fname in data_dict:
    data_dict[fname] = load_file(fname)

def plot(data: list, fname: str):
  plt.plot(data[0], data[1])
  plt.title("Integral of " + fname)
  plt.savefig("/gluster/data/atlas/jdombrowski/SelectionCutIntegrals/" + fname + ".png")

def get_sig(data_dict: dict, signal: str, background: str):
  # signal and background must be same size
  n = len(data_dict[signal][0])
  if n != len(data_dict[background][0]):
    print("Signal and Bckgrnd not of same size", n, len(data_dict[background][0]))
    return -1
  sig = [data_dict[signal][0], [0 for i in range(0, n)]]
  for i in range(0, n):
    if(data_dict[background][1][i] == 0): 
      sig[1][i] = 0
      continue
    sig[1][i] = data_dict[signal][1][i]/np.sqrt(data_dict[background][1][i])
  data_dict[signal + " significance"] = sig
  return 0

def plot_all(data_dict):
  for fname in data_dict:
    print("Loading...",fname)
    plot(data_dict[fname], fname)
    plt.show()


def main():
  data_dict = {
    "DeltaPhiIn DeltaRIndicies METbetweenJets DeltaRCut0p15" : [[], []],
    "DeltaPhiOut DeltaRIndicies METbetweenJets DeltaRCut0p15": [[], []]
  }
  load_dict(data_dict)

  get_sig(data_dict, "DeltaPhiIn DeltaRIndicies METbetweenJets DeltaRCut0p15", "DeltaPhiOut DeltaRIndicies METbetweenJets DeltaRCut0p15")
  plot_all(data_dict)
  # plot(data_dict["DeltaPhiIn DeltaRIndicies METbetweenJets DeltaRCut0p15 significance"], "sig")
  return 0;

if __name__ == "__main__":
  main()