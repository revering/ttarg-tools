 #!/usr/bin/env python

import argparse
from time import strftime
import os
import math
import ROOT
from array import array

parser = argparse.ArgumentParser(description="Parse Root Files to create angle and energy distributions")
parser.add_argument("rootfile", help="name of files in directory")
parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd())
arg = parser.parse_args()

outDir = arg.outDir
nevts = 0
#Check for trailing slash on ouput dir and delete
if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
infileType = "."+arg.rootfile.split(".")[-1]

if not os.path.isdir(outDir):
     print "Output directory does not exist!"
     quit()

if not os.path.isfile(arg.rootfile):
     print "Input file does not exist!"
     quit()

if infileType != ".root":
     print "Input file is of incorrect type \"%s\"!"%(infileType)
     quit()

rfile = ROOT.TFile(arg.rootfile,"UPDATE")
tree = rfile.Get("vector_tree")
pvec = ROOT.TLorentzVector()
mvec = ROOT.TLorentzVector()

tree.SetBranchAddress("mu_pos_vec",pvec)
tree.SetBranchAddress("mu_neg_vec",mvec)

entries = tree.GetEntries()
cut_direc = rfile.mkdir("cut_events")
cut_direc.cd()
otree = ROOT.TTree("cut_events","Tree containing lorentz vectors after cuts")
otree.Branch("mu_pos_vec","TLorentzVector",pvec)
otree.Branch("mu_neg_vec","TLorentzVector",mvec)

for j in range(entries):
   tree.GetEntry(j)
   if abs(pvec.Eta())<2.4 and abs(mvec.Eta())<2.4:
      if pvec.Pt()>20 and mvec.Pt()>20:
         invmass = pvec + mvec
	 if invmass.M()>60 and invmass.M()<120:
	    otree.Fill()

otree.Write()
rfile.Close()
