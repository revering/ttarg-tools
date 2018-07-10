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

# inputs: root file
  # all_events TTree
    # mu_pos
    # mu_neg
# ouputs: root file with histograms
  # all_events TTree
    # mu_pos
    # mu_neg
  # histograms
    # pt
    # eta
    # phi
    # inv_m
file1 = sys.argv[1]

class Plotter:
  # Class which takes in a root file and plots pt, eta, phi, and p^2 of the particles, designed for muons from Z->dimuon events
  def __init__(self,TFile_in):
    self.tfin = TFile_in

def plot(f1):
  # check file
  fin1 = ROOT.TFile(f1,"UPDATE")
  if not fin1.IsOpen():
    print "Could not open {0}".format(f1)
    return
  # access tree
  all_tree = fin1.Get("all_events")
  has_histograms = fin1.GetDirectory("histograms")
  if has_histograms:
    print "found histograms"
  else:
    print "didn't find histograms, making histograms directory..."
    histograms = fin1.mkdir("histograms")
  fin1.cd("histograms")
  mu_pos = ROOT.TLorentzVector()
  mu_neg = ROOT.TLorentzVector()
  all_tree.SetBranchAddress("mu_pos_vec",mu_pos)
  all_tree.SetBranchAddress("mu_neg_vec",mu_neg)
  ptpos = ROOT.TH1F("ptpos","Transverse Momentum (mu+)",200,0,100)
  ptneg = ROOT.TH1F("ptneg","Transverse Momentum (mu-)",200,0,100)
  etapos = ROOT.TH1F("etapos","#eta (mu+)",200,-100,100)
  etaneg = ROOT.TH1F("etaneg","#eta (mu-)",200,-100,100)
  phipos = ROOT.TH1F("phipos","#phi (mu+)",200,-100,100)
  phineg = ROOT.TH1F("phineg","#phi (mu-)",200,-100,100)
  invm = ROOT.TH1F("invm","Invariant Mass (mu+)",200,-100,100)
  for i in range(all_tree.GetEntries()):
    all_tree.GetEntry(i)
    ptpos.Fill(mu_pos.Pt())
    ptneg.Fill(mu_neg.Pt())
    etapos.Fill(mu_pos.Eta())
    etaneg.Fill(mu_neg.Eta())
    phipos.Fill(mu_pos.Phi())
    phineg.Fill(mu_neg.Phi())
    invm.Fill(mu_pos.M())
    invm.Fill(mu_neg.M())
  #ROOT.gStyle.SetOptStat(0)
  canvas = ROOT.TCanvas("canvas","Z->mumu",1400,700)
  canvas.Divide(4,2)
  canvas.cd(1)
  ptpos.Draw()
  canvas.cd(2)
  ptneg.Draw()
  canvas.cd(3)
  etapos.Draw()
  canvas.cd(4)
  etaneg.Draw()
  canvas.cd(5)
  phipos.Draw()
  canvas.cd(6)
  phineg.Draw()
  canvas.cd(7)
  invm.Draw()
  fin1.Write()
  #name1 = pytools.find_jobname(f1)
  #name2 = pytools.find_jobname(f2)
  # legend
  #leg = ROOT.TLegend(0.6,0.7,0.9,0.9)
  #leg.AddEntry(hist1,name1)
  #leg.AddEntry(hist2,name2)
  #leg.Draw()
  raw_input("waiting...")
  return
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
plot(file1)
