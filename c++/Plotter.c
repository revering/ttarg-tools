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
}

void Plotter::book(TFile* save_dir, const std::string sub_dir) {
  save_dir->mkdir(sub_dir.c_str());
  save_dir->cd(sub_dir.c_str());
  pt_pos = new TH1F("pt_pos","Transverse Momentum (mu+)",100,0,100);
  pt_neg = new TH1F("pt_neg","Transverse Momentum (mu-)",100,0,100);
  eta_pos = new TH1F("eta_pos","#eta (mu+)",160,-10,10);
  eta_neg = new TH1F("eta_neg","#eta (mu-)",160,-10,10);
  phi_pos = new TH1F("phi_pos","#phi (mu+)",100,-TMath::Pi(),TMath::Pi());
  phi_neg = new TH1F("phi_neg","#phi (mu+)",100,-TMath::Pi(),TMath::Pi());
  inv_m = new TH1F("inv_m","Invariant Mass",150,0,150);
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
}

void Plotter::write(TDirectory* save_dir) {
  save_dir->Write();
}
