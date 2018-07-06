#!/usr/bin/env python

# author: Reese Petersen, University of Minnesota
# email: pet00831@umn.edu

import os
import sys
import argparse
import ROOT

# assumptions: the lhe file follows the Les Houches Event File standard, with a row of 13 numbers for each particle in each event, 
# the first being the pdgID and the seventh - tenth being the px, py, pz, and energy, respectively

dir_out_d = "/home/%s/"%(os.environ["USER"])
parser = argparse.ArgumentParser(description = "Extracts energy-momentum 4-vectors of mu+ and mu- in lhe files and saves them in a Ttree in a root file")
parser.add_argument("-fi","--file_in", dest = "file_in", help = "lhe file to extract di-lepton 4-vectors from", required = True)
parser.add_argument("-do","--directory_out", dest = "directory_out", help = "lhe file to extract di-lepton 4-vectors from", default = dir_out_d)
arg = parser.parse_args()
# check the input file
flhe = arg.file_in
if not os.path.isfile(flhe):
  print "Zll_lhe_to_Ttrees.py: %s not found. Exiting..."%(flhe)
  quit()
if flhe.split('.')[-1] != "lhe":
  print "Zll_lhe_to_Ttrees.py: %s does not have a '.lhe' extension."%(flhe)
  choice = raw_input("Continue anyway? (y/n)")
  if choice != "y":
    print "Zll_lhe_to_Ttrees.py: Exiting..."
    quit()
# make the output directory if it does not exist
if arg.directory_out.split('/')[-1] != "":
  out_dir = arg.directory_out+"/"
else:
  out_dir = arg.directory_out
if os.path.exists(out_dir):
  print "Zll_lhe_to_Ttrees.py: %s exists"%(out_dir)
else:
  print "Zll_lhe_to_Ttrees.py: %s does not exist, creating..."%(out_dir)
  os.makedirs(out_dir)
# initialize a Ttree
if flhe.find('/') > -1:
  lhename = flhe.split('/')[-1]
else:
  lhename = flhe
outname = lhename.split('.lhe')[0]
fout = ROOT.TFile(out_dir+"Zll_4vec_"+outname+".root")
ROOT.Ttree()
# read out the energy and momenta from the lhe file
lhe = open(flhe,"r")
curline = ""
while curline.find("<event>") == -1:
  curline = lhe.readline()
the_end = False
end_loop_count = 0
event_loop = 0
muon_count = 0
for line in lhe:
  curlist = line.split(' ')
  for i in range(curlist.count('')):
    curlist.remove('')
  if len(curlist) == 13 and abs(int(curlist[0])) == 13:
    muon_count += 1
    #pxs = curlist[6]
    #pys = curlist[7]
    #pzs = curlist[8]
    #ens = curlist[9]
print "Zll_lhe_to_Ttrees.py: total number of mu+/-: %d"%(muon_count)
quit()
