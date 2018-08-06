#include "Plotter.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TFile.h"
#include <string>
#include "Event.h"
#include "TMath.h"

void Plotter::book(TDirectory* save_dir, const std::string sub_dir) {
  save_dir->mkdir(sub_dir.c_str());
  save_dir->cd(sub_dir.c_str());
  pt_pos = new TH1F("pt_pos","Transverse Momentum (mu+)",100,0,100);
  pt_neg = new TH1F("pt_neg","Transverse Momentum (mu-)",100,0,100);
  eta_pos = new TH1F("eta_pos","#eta (mu+)",160,-10,10);
  eta_neg = new TH1F("eta_neg","#eta (mu-)",160,-10,10);
  phi_pos = new TH1F("phi_pos","#phi (mu+)",100,-TMath::Pi(),TMath::Pi());
  phi_neg = new TH1F("phi_neg","#phi (mu+)",100,-TMath::Pi(),TMath::Pi());
  inv_m = new TH1F("inv_m","Invariant Mass",150,0,150);
  pt_change = new TH1F("pt_change","Change in pt",250,0,200);
  inv_pt_change = new TH1F("inverse_pt_change","Inverse Pt Change",200,0,1);
  pt_brem = new TH1F("pt_brem","Transverse Momentum Of Brem Muon",100,0,100);
  pt_ratio = new TH1F("pt_ratio","Brem Pt Ratio",200,0,1);
  eta_brem = new TH1F("eta_brem","#eta (mu+)",160,-10,10);
  int_pt_ratio = new TH1F("int_pt_ratio", "Integrated Pt Ratio", 200,0,1);
}


void Plotter::book(TFile* save_dir, const std::string sub_dir) {
  save_dir->mkdir(sub_dir.c_str());
  save_dir->cd(sub_dir.c_str());
  root_dir = save_dir;
  hist_dir = sub_dir;
  pt_pos = new TH1F("pt_pos","Transverse Momentum (mu+)",100,0,100);
  pt_neg = new TH1F("pt_neg","Transverse Momentum (mu-)",100,0,100);
  eta_pos = new TH1F("eta_pos","#eta (mu+)",160,-10,10);
  eta_neg = new TH1F("eta_neg","#eta (mu-)",160,-10,10);
  phi_pos = new TH1F("phi_pos","#phi (mu+)",100,-TMath::Pi(),TMath::Pi());
  phi_neg = new TH1F("phi_neg","#phi (mu+)",100,-TMath::Pi(),TMath::Pi());
  inv_m = new TH1F("inv_m","Invariant Mass",150,0,150);
  pt_change = new TH1F("pt_change","Change in pt", 250,0,200);
  pt_brem = new TH1F("pt_brem","Transverse Momentum of Brem Muon",100,0,100);
  inv_pt_change = new TH1F("inverse_pt_change","Inverse Pt Change",200,0,1);
  pt_ratio = new TH1F("pt_ratio","Brem Pt Ratio",200,0,1);
  eta_brem = new TH1F("eta_brem","#eta (mu+)",160,-10,10);
  int_pt_ratio = new TH1F("int_pt_ratio", "Integrated Pt Ratio", 200,0,1);
}

void Plotter::fill(event e) {
  pt_pos->Fill(e.mu_pos_vec.Pt());
  pt_neg->Fill(e.mu_neg_vec.Pt());
  eta_pos->Fill(e.mu_pos_vec.Eta());
  eta_neg->Fill(e.mu_neg_vec.Eta());
  phi_pos->Fill(e.mu_pos_vec.Phi());
  phi_neg->Fill(e.mu_neg_vec.Phi());
  TLorentzVector sum = e.mu_pos_vec+e.mu_neg_vec;
  inv_m->Fill(sum.M());
  if(e.pbrem==1)
  {
     pt_brem->Fill(e.mu_pos_vec.Pt());
     eta_brem->Fill(e.mu_pos_vec.Eta());
     pt_change->Fill(e.pre_brem_vec.Pt()-e.mu_pos_vec.Pt());
     inv_pt_change->Fill(-1./e.pre_brem_vec.Pt()+1./e.mu_pos_vec.Pt());
     pt_ratio->Fill(e.mu_pos_vec.Pt()/e.pre_brem_vec.Pt());
  }
  else if(e.pbrem==-1)
  {
     pt_brem->Fill(e.mu_neg_vec.Pt());
     eta_brem->Fill(e.mu_neg_vec.Eta());
     pt_change->Fill(e.pre_brem_vec.Pt()-e.mu_neg_vec.Pt());
     inv_pt_change->Fill(1./e.mu_neg_vec.Pt()-1./e.pre_brem_vec.Pt());
     pt_ratio->Fill(e.mu_neg_vec.Pt()/e.pre_brem_vec.Pt());
  }
}

void Plotter::integrate(TH1* inthist, TH1* initial)
{
   Int_t nbins = initial->GetNbinsX();
   Int_t sum = 0;
   for(int i=0;i<=nbins;i++)
   {
      inthist->SetBinContent(i,initial->GetBinContent(i)+sum);
      sum = inthist->GetBinContent(i);
   }
   if(initial->GetEntries()!=0) 
   {
      inthist->Scale(1/initial->GetEntries());
      inthist->SetEntries(initial->GetEntries());
   }
   return;
}

void Plotter::write() 
{
//   root_dir->cd(hist_dir.c_str());
   integrate(int_pt_ratio, pt_ratio);
}
