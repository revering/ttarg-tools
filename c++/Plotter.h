#ifndef plotter_h
#define plotter_h
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include <string>
#include "Event.h"

class Plotter {
  public:
    void book(TDirectory* save_dir, const std::string sub_dir);
    void book(TFile* save_dir, const std::string sub_dir);
    TDirectory* root_dir;
    std::string hist_dir;
    void integrate(TH1* inthist, TH1* initial);
    void fill(event e);
    void write();
    void fill_int();
    TH1* pt_pos;
    TH1* pt_neg;
    TH1* eta_pos;
    TH1* eta_neg;
    TH1* phi_pos;
    TH1* phi_neg;
    TH1* inv_m;
    TH1* pt_change;
    TH1* pt_brem;
    TH1* eta_brem;
    TH1* inv_pt_change;
    TH1* pt_ratio;
    TH1* int_pt_ratio;
};

#endif
