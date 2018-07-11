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

#endif
