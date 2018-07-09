#!/usr/bin/env python

# import useful libraries
import os
import math
import sys
import argparse
import subprocess
import ROOT
import datetime
import pytools

# taking an argument from the command line
#parser = argparse.ArgumentParser(description = 'This plots a monte carlo simulation for the energy flux of cosmic muons at sea level on Earth as a function of energy.')
#parser.add_argument("--rootFiles" , dest="rootFiles" , help="The root files containing histograms you want to plot")
#arg = parser.parse_args()

file1 = sys.argv[1]
#file2 = sys.argv[2]

def plot(f1):#,f2):
  fin1 = ROOT.TFile(f1)
  #fin2 = ROOT.TFile(f2)
  if not fin1.IsOpen(): #or not fin2.IsOpen():
    print "Could not open" #one or both files, exiting..."
    return
  # histograms
  tree = fin1.Get("vector_tree")
  mu_pos = ROOT.TLorentzVector()
  mu_neg = ROOT.TLorentzVector()
  tree.SetBranchAddress("mu_pos_vec",mu_pos)
  tree.SetBranchAddress("mu_neg_vec",mu_neg)
  pxpos = ROOT.TH1F("pxp","pxp",200,-100,100)
  pypos = ROOT.TH1F("pyp","pyp",200,-100,100)
  pzpos = ROOT.TH1F("pzp","pzp",200,-100,100)
  pxneg = ROOT.TH1F("pxn","pxn",200,-100,100)
  pyneg = ROOT.TH1F("pyn","pyn",200,-100,100)
  pzneg = ROOT.TH1F("pzn","pzn",200,-100,100)
  for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    pxpos.Fill(mu_pos.Px())
    pypos.Fill(mu_pos.Py())
    pzpos.Fill(mu_pos.Pz())
    pxneg.Fill(mu_neg.Px())
    pyneg.Fill(mu_neg.Py())
    pzneg.Fill(mu_neg.Pz())
  ROOT.gStyle.SetOptStat(0)
  canvas = ROOT.TCanvas("canvas","Z->mumu",900,700)
  canvas.Divide(3,2)
  canvas.cd(1)
  pxpos.Draw()
  canvas.cd(2)
  pypos.Draw()
  canvas.cd(3)
  pzpos.Draw()
  canvas.cd(4)
  pxneg.Draw()
  canvas.cd(5)
  pyneg.Draw()
  canvas.cd(6)
  pzneg.Draw()
  #hist2 = fin2.Get("analyzer/Energy")
  #hist1.SetLineColor(1)
  #hist2.SetLineColor(2)
  #hist1.Draw()
  #hist2.Draw("SAME")
  # job name
  #name1 = pytools.find_jobname(f1)
  #name2 = pytools.find_jobname(f2)
  # legend
  #leg = ROOT.TLegend(0.6,0.7,0.9,0.9)
  #leg.AddEntry(hist1,name1)
  #leg.AddEntry(hist2,name2)
  #leg.Draw()
  raw_input("waiting...")
  return
  #fin.ls()

  #ROOT.TIter current(fin.GetListOfKeys())
  
  #fin.GetListOfKeys().Print()
  #ROOT.gDirectory.GetListOfKeys()
  begin = datetime.datetime.now()
  print "start: "+str(begin.hour)+":"+str(begin.minute)+":"+str(begin.second)
  # checking for the output directory
  if output != "":
    if not os.path.exists(output):
      print "Output directory does not exist! Exiting..."
      quit()
  
  # make sure there is a slash at the end of the output path
  if output.split('/')[-1] != '':
    outputpath = output+'/'
  # canvas
  pycan = ROOT.TCanvas("canvy","I am a canvas title",1700,900)
  pycan.SetLogx();
  # legend
  pylegend = ROOT.TLegend(0.7,0.7,0.9,0.9)
plot(file1)#,file2)
