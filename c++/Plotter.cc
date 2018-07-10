#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include <string>
#include "Event.h"

class Plotter {
  public:
    void book(TDirectory* save_dir, const string sub_dir);
    void fill(event e);
    void write(TDirectory* save_dir);
    TH1* pt_pos;
    TH1* pt_neg;
    TH1* eta_pos;
    TH1* eta_neg;
    TH1* phi_pos;
    TH1* phi_neg;
    TH1* inv_m;
};

void Plotter::book(TDirectory* save_dir, const string sub_dir) {
  save_dir->mkdir(sub_dir.c_str());
  save_dir->cd(sub_dir.c_str());
  pt_pos = new TH1F("pt_pos","Transverse Momentum (mu+)",100,0,100);
  pt_neg = new TH1F("pt_neg","Transverse Momentum (mu-)",100,0,100);
  eta_pos = new TH1F("eta_pos","#eta (mu+)",100,-50,50);
  eta_neg = new TH1F("eta_neg","#eta (mu-)",100,-50,50);
  phi_pos = new TH1F("phi_pos","#phi (mu+)",360,0,360);
  phi_neg = new TH1F("phi_neg","#phi (mu+)",360,0,360);
  inv_m = new TH1F("inv_m","Invariant Mass",200,0,200);
}

void Plotter::fill(event e) {
  pt_pos->Fill(e.mu_pos_vec.Pt());
  pt_neg->Fill(e.mu_neg_vec.Pt());
  eta_pos->Fill(e.mu_pos_vec.Eta());
  eta_neg->Fill(e.mu_neg_vec.Eta());
  phi_pos->Fill(e.mu_pos_vec.Phi());
  phi_neg->Fill(e.mu_neg_vec.Phi());
  inv_m->Fill(e.mu_pos_vec.M());
  inv_m->Fill(e.mu_neg_vec.M());
}

void Plotter::write(TDirectory* save_dir) {
  save_dir->Write();
}
