import ROOT as RT
from glob import glob as ls
from argparse import ArgumentParser as ap

parser = ap()

parser.add_argument( "-p", type=str, help='Path to files', default='./')
args = parser.parse_args()

#print(ls(args.p))

all_files = [i.split("/")[-1] for i in ls(args.p + "*")]
files = []
global_meds = []
for p in all_files:
  if 'global_median' in p:
    print(p)
    global_meds.append(p)
runs = set([i.replace(".txt", "").split("_")[-1] for i in global_meds])
print(runs)

for r in runs:
  f0 = open(args.p + 'global_median_0_' + r + '.txt', 'r') 
  f1 = open(args.p + 'global_median_1_' + r + '.txt', 'r') 
  f2 = open(args.p + 'global_median_2_' + r + '.txt', 'r') 

  info_0 = f0.readlines()[0]
  info_1 = f1.readlines()[0]
  info_2 = f2.readlines()[0]
  print(info_0)
  print(info_1)
  print(info_2)

  run_num = info_0.split()[0]
  norm_0 = info_0.split()[1]
  norm_err_0 = str(float(norm_0)/10.)
  print(run_num, norm_0, norm_err_0)

  norm_1 = info_1.split()[1]
  norm_err_1 = str(float(norm_1)/10.)

  norm_2 = info_2.split()[1]
  norm_err_2 = str(float(norm_2)/10.)


  output = [
    'channel,tv,norm,norm_err\n',
    ','.join(['0', run_num, norm_0, norm_err_0+"\n" ]),
    ','.join(['1', run_num, norm_1, norm_err_1+"\n" ]),
    ','.join(['2', run_num, norm_2, norm_err_2+"\n" ])
  ]
  print(output)
  outfile = open('globalmedians_run' + run_num + '.csv', 'w')
  outfile.writelines(output)
  outfile.close()
