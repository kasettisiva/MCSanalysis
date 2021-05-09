import os
from argparse import ArgumentParser as ap

parser = ap()

parser.add_argument( "-r", type=str, help='Run number', default='5387')
#parser.add_argument( "-o", type=str, help='Time (s) to sleep between calls', default=)

args = parser.parse_args()

fout = open('file_list_' + args.r + '.txt', 'w')


getnames = os.popen("samweb list-definition-files protodune-sp_runset_" + args.r + "_michelremoving_merged_v09_09_01_v0")
filenames = getnames.readlines()
out_names = []
for fn in filenames:
  fn = fn.rstrip()
  fileloc = os.popen(f"samweb locate-file {fn}")
  fileloc = fileloc.readlines()
  fileloc = fileloc[0].split(":")
  fileloc = fileloc[1].split("(")
  print (fileloc[0]+"/"+fn)
  out_names.append(fileloc[0]+"/"+fn +"\n")

fout.writelines(out_names)
