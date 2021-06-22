import ROOT as RT
from argparse import ArgumentParser as ap
from array import array

parser = ap()

parser.add_argument("-i", type=str, help='Input file', default = "")
parser.add_argument("-o", type=str, help='Output file', default = "covariance.root")
args = parser.parse_args()

input_file = open(args.i, "r")

output_file = RT.TFile(args.o, "RECREATE")

lines = input_file.readlines()
cov_vals = []
n = 0
for l in lines:
  if "#" in l:
    print(l)
    continue 
  l = l.strip("\n")
  l_split = l.split()
  print(l_split)
  #cov_vals.append([])
  n += 1
  for s in l_split:
#    cov_vals[-1].append(float(s))
    cov_vals.append(float(s))

cov_matrix = RT.TMatrixD(n, n, array("d", cov_vals))
cov_matrix.Write("m")
output_file.Close()
