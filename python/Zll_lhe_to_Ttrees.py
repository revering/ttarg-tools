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
parser.add_argument("-fi","--files_in", dest = "files_in", help = "lhe file(s) to extract di-lepton 4-vectors from", required = True, nargs="+")
parser.add_argument("-do","--directory_out", dest = "directory_out", help = "output directory, default = %s"%(dir_out_d), default = dir_out_d)
arg = parser.parse_args()
# check the input file
flhe = arg.files_in
flhe0 = flhe[0]
if not os.path.isfile(flhe0):
  print "Zll_lhe_to_Ttrees.py: %s not found. Exiting..."%(flhe0)
  quit()
if flhe0.split('.')[-1] != "lhe":
  print "Zll_lhe_to_Ttrees.py: %s does not have a '.lhe' extension."%(flhe0)
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
# open a TFile, initialize a TTree and TLorentzVectors
if flhe0.find('/') > -1:
  lhename = flhe0.split('/')[-1]
else:
  lhename = flhe0
outname = "_".join((lhename.split('.lhe')[0]).split('_')[:-1])
fout = ROOT.TFile(out_dir+"Zll_4vec_"+outname+".root","CREATE")
tree = ROOT.TTree("all_events", "Tree of mu+- 4-vectors (Z->mu+-)")
mu_pos_vec = ROOT.TLorentzVector()
mu_neg_vec = ROOT.TLorentzVector()
tree.Branch("mu_pos_vec","TLorentzVector",mu_pos_vec)
tree.Branch("mu_neg_vec","TLorentzVector",mu_neg_vec)
# read out the energy and momenta from the lhe file
for infile in flhe:
  lhe = open(infile,"r")
  for line in lhe:
    curlist = line.split(' ')
    for i in range(curlist.count('')):
      curlist.remove('')
    if len(curlist) == 13 and (curlist[0] == "13" or curlist[0] == "-13"):
      pxs = curlist[6]
      pys = curlist[7]
      pzs = curlist[8]
      ens = curlist[9]
      if int(curlist[0]) == 13: #mu+
        mu_pos_vec.SetPxPyPzE(float(pxs),float(pys),float(pzs),float(ens))
      elif int(curlist[0]) == -13: #mu-
        mu_neg_vec.SetPxPyPzE(float(pxs),float(pys),float(pzs),float(ens))
    if line.find("</event>") > -1:
      tree.Fill()
  lhe.close()
#fout.Write()
fout.Close()
quit()
