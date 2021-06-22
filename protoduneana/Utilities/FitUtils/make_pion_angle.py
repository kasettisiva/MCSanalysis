import ROOT as RT
from argparse import ArgumentParser as ap
from math import exp
parser = ap()

parser.add_argument("-i", type=str, help='Input file', default = "")
parser.add_argument("-o", type=str, help='Output file', default = "angle_weight.root")
parser.add_argument("-t", type=int, help='Shift type', default=0)
parser.add_argument("-p", type=float, help='Position -- only for gaus', default=.5)
args = parser.parse_args()

f = RT.TFile(args.i, "OPEN")
t = f.Get("pduneana/beamana")

end_p = [0., .4, .6,   .8,   1.]
nbins = [  5,   25, 25,  25,  10]
shifts = [ 2,  12, 12, 12, 5]

RT.gROOT.SetBatch() 
fout = RT.TFile(args.o, "RECREATE")

rs = []

def shift_bins(n_bins, shift):
  bins = range(1, n_bins+1)
  new_bins = [i + shift for i in range(1, len(bins)-shift+1)] + [i - len(bins) + shift for i in range(len(bins)-shift+1, len(bins)+1)]
  print(new_bins)
  return new_bins
for i in range(1, 12):
  shift_bins(10, i)

def flip_pions(h_in):
  h2 = h_in.Clone()
  for j in range(1, h2.GetNbinsX()+1):
    h2.SetBinContent(j, h.GetBinContent(h2.GetNbinsX() + 1 - j))
  return h2

def flat(h_in):
  h2 = h_in.Clone()
  for j in range(1, h2.GetNbinsX()+1): 
    h2.SetBinContent(j, 1.)
  h2.Scale(1. / h2.Integral())
  return h2

def shift_pions(h_in, i):
  new_bins = shift_bins(h_in.GetNbinsX(), shifts[i])
  h2 = h_in.Clone()
  for j in range(1, h2.GetNbinsX()+1):
    h2.SetBinContent(j, h.GetBinContent(new_bins[j-1]))
  h2.Scale(1. / h2.Integral())
  return h2

def gaus_pions(h_in):
  n_bins = h_in.GetNbinsX()
  width = n_bins/4.
  mean = n_bins*args.p
  h2 = h_in.Clone()
  for i in range(1, n_bins+1):
    h2.SetBinContent(i, exp(-.5*(((i - mean)/width)**2)))
  h2.Scale(1. / h2.Integral())
  return h2 

def check_bins(h_in):
  for i in range(1, h_in.GetNbinsX()+1):
    if h_in.GetBinContent(i) == 0.:
      print(h_in.GetName(), "empty:", i)

def get_variation(h_in, t, j = 0):
  if t == 0:
    return flip_pions(h_in)
  elif t == 1:
    return flat(h_in)
  elif t == 2:
    return shift_pions(h_in, j)
  elif t == 3:
    return gaus_pions(h_in)
  else:
    print("error, unkown type", t)

for i in range(1, len(end_p)):
  print("hist", i, "nbins:", nbins[i-1])
  h = RT.TH1D("h" + str(i), "", nbins[i-1], -1., 1.)
  t.Draw("(true_beam_daughter_startPz*true_beam_endPz + true_beam_daughter_startPx*true_beam_endPx + true_beam_daughter_startPy*true_beam_endPy)/(true_beam_daughter_startP*true_beam_endP)" + 
         ">>h" + str(i),
         #">>h" + str(i) + "(" + str(nbins[i-1]) + "-1., 1.)",
         "true_daughter_nPiPlus == 1 && true_daughter_nPiMinus == 0 && " +
         "new_interaction_topology == 3 && true_beam_PDG == 211 && abs(true_beam_daughter_PDG) == 211 && true_beam_endP >= " + str(end_p[i-1]) + " && true_beam_endP < " + str(end_p[i]))
  #h = RT.gDirectory.Get("h" + str(i))
  h.Scale(1./h.Integral())
  h.Write()

  #if args.t == 0:
  #  h2 = flip_pions(h)
  #else:
  #  h2 = flat(h)
  h2 = get_variation(h, args.t, i-1)
  print(h2.Integral())
  r = h2.Clone()
  r.Divide(h)
  h2.Write("h" + str(i) + "_rev")
  r.Write("r" + str(i))
  rs.append(r)
  
h = RT.TH1D("h" + str(len(end_p)), "", nbins[-1], -1., 1.)
t.Draw("(true_beam_daughter_startPz*true_beam_endPz + true_beam_daughter_startPx*true_beam_endPx + true_beam_daughter_startPy*true_beam_endPy)/(true_beam_daughter_startP*true_beam_endP)" + 
       ">>h" + str(len(end_p)),
       #">>h" + str(len(end_p)) + "(" + str(nbins[-1]) + "-1., 1.)",
       "true_daughter_nPiPlus == 1 && true_daughter_nPiMinus == 0 && " +
       "new_interaction_topology == 3 && true_beam_PDG == 211 && abs(true_beam_daughter_PDG) == 211 && true_beam_endP >= " + str(end_p[-1]))
#h = RT.gDirectory.Get("h" + str(len(end_p)))
h.Scale(1./h.Integral())
h.Write()
#if args.t == 0:
#  h2 = flip_pions(h)
#else:
#  h2 = flat(h)
h2 = get_variation(h, args.t, len(end_p)-1)
print(h2.Integral())
r = h2.Clone()
r.Divide(h)
h2.Write("h" + str(len(end_p)) + "_rev")
r.Write("r" + str(len(end_p)))
rs.append(r)


fout.Close()

