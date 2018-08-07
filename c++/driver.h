/**
 * @file driver.h
 * @brief Class driving the dark brem simulation + analysis.
 * @author Michael Revering, University of Minnesota
 */

#ifndef driver_h
#define driver_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "DarkPhotons.hh"
#include "Plotter.h"
// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"
#include "Event.h"

//Class to load in tfiles, store histograms, and hold functions for performing analysis.
class driver {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   double          map; //A' mass
   DarkPhotons*    dphoton;
   Plotter         all_events;
   Plotter         inv_mass_cut;
   Plotter         after_cuts;
   Plotter         after_brem;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   TLorentzVector  *mu_pos_vec;
   TLorentzVector  *mu_neg_vec;

   // List of branches
   TBranch        *b_mu_pos_vec;   //!
   TBranch        *b_mu_neg_vec;   //!

   driver(TTree *tree=0, double mAPrime=1.0, std::string ifname = "");
   virtual ~driver();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   int              Select_Lepton(event Evnt);
   TLorentzVector   simulate_dbrem(TLorentzVector initial);
   bool             pass_cuts(event evnt);
   bool             eta_cuts(TLorentzVector muon, TLorentzVector track);
   bool             pt_cuts(TLorentzVector muon, TLorentzVector track);
   bool             dot_cuts(TLorentzVector muon, TLorentzVector track);
   bool             pass_inv_mass(event evnt);
   event            GetEvent(double entry);
   double           GetEntries();
};

#endif

#ifdef driver_cxx
driver::driver(TTree *tree, double mAPrime, std::string ifname) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Zll_4vec_cmsgrid_final_RandomSeed_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Zll_4vec_cmsgrid_final_RandomSeed_0.root");
      }
      f->GetObject("all_events",tree);

   }
   Init(tree);
   map = mAPrime;
   dphoton = new DarkPhotons(map, 0., 1., 28.,14., 2.32, 1., ifname);
   dphoton->PrepareTable();
}

driver::~driver()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t driver::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t driver::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void driver::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mu_pos_vec = 0;
   mu_neg_vec = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mu_pos_vec", &mu_pos_vec, &b_mu_pos_vec);
   fChain->SetBranchAddress("mu_neg_vec", &mu_neg_vec, &b_mu_neg_vec);
   Notify();
}

Bool_t driver::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void driver::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t driver::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

int driver::Select_Lepton(event evnt)
{
// Selects which of the two leptons undergoes a dark brem. Weights the choice by the cross
// section of each particle, which differs due to their energies. Returns a 1 if the positive
// lepton brems, and a -1 if the negative one does. Returns 0 if both cross sections are zero.
  double MuPlusE = evnt.mu_pos_vec.E();
  double MuMinusE = evnt.mu_neg_vec.E();
  double pWeight = dphoton->GetsigmaTot(MuPlusE);
  double mWeight = dphoton->GetsigmaTot(MuMinusE);
  if ((pWeight+mWeight)==0)
  {
     //Can happen if the energies are too low, or out of range of the table.
     std::cerr << "Cross section is zero for both leptons.\n"; 
     return 0;
  }
  if((pWeight<0))
  {
     //Happens occasionally due to a bad function interpolation. 
     printf("Negative cross section, energy out of range, Failing E is %e\n",MuPlusE);
     return 0;
  }
  if((mWeight<0))
  { 
     printf("Negative cross section, energy out of range, Failing E is %e\n",MuMinusE);
     return 0;
  }

  //Make the selection.
  double A = drand48()*(pWeight+mWeight);
  if (A<=pWeight) {return 1;}
  else {return -1;}
}

TLorentzVector   driver::simulate_dbrem(TLorentzVector initial)
{
   //Simulates the dark brem. Takes in the 4-vector of the incoming muon, returns 
   //the 4-vector of the outgoing muon after the brem.
   double energy = initial.E();
   TVector3 UzVec = initial.Vect().Unit();
   TLorentzVector* outvec = (TLorentzVector*)dphoton->MuSimulateEmission(energy)->Clone();
   outvec->RotateUz(UzVec);
   return *outvec;
}

bool driver::pass_inv_mass(event evnt)
{
//Invariant mass cut. Done before choosing the brem to avoid save simulation time and to
//avoid events with zero cross sections.
   bool pass = true;
   TLorentzVector invmass = evnt.mu_neg_vec+evnt.mu_pos_vec;
   if ((invmass.M()<81.2)||(invmass.M()>101.2)) {pass=false;}
   return pass;
}

bool driver::pass_cuts(event evnt)
{
//Make the cuts for the muon+track. Set to match those made for the disappearing tracks search.
//"muon" refers the the lepton which does not brem, "track" refers to the one that does.
   bool pass = true;
   TLorentzVector muon, track;
   if(evnt.pbrem==1)
   {
      muon = evnt.mu_neg_vec;
      track = evnt.pre_brem_vec;
   }
   else if(evnt.pbrem==-1)
   {
      muon = evnt.mu_pos_vec;
      track = evnt.pre_brem_vec;
   }
   else {pass=false;}
   if(pt_cuts(muon,track)==false) {pass=false;}
   else if(eta_cuts(muon,track)==false) {pass=false;}
   else if(dot_cuts(muon,track)==false) {pass=false;} 
   return pass;
}

bool driver::pt_cuts(TLorentzVector muon, TLorentzVector track)
{
   bool pass = true;
   if(muon.Pt()<25) {pass=false;}
   else if(track.Pt()<30) {pass=false;}
   return pass;
}

bool driver::eta_cuts(TLorentzVector muon, TLorentzVector track)
{
//Eta cuts. The abs() function wasn't behaving as I expected, so needed to find the 
//absolute value manually.
   bool pass = true;
   double m_eta = muon.Eta();
   double t_eta = track.Eta();
   if (t_eta<0) {t_eta=t_eta*-1.;}
   if (m_eta<0) {m_eta=m_eta*-1.;}
   if(m_eta>2.1) {pass=false;}
   else if(t_eta>2.1) {pass=false;}
   else if((t_eta>0.15)&&(t_eta<0.35)) {pass=false;}
   else if((t_eta>1.42)&&(t_eta<1.65)) {pass=false;}
   else if((t_eta>1.55)&&(t_eta<1.85)) {pass=false;}
   return pass;
}

bool driver::dot_cuts(TLorentzVector muon, TLorentzVector track)
{
//Cut on the dot product of the two 3-momenta. Makes sure the angle between the
//muon and the track is at least 90 degrees.
   bool pass = true;
   TLorentzVector z_bos = muon+track;
   double p1dotp2 = muon.X()*track.X()+muon.Y()*track.Y()+muon.Z()*track.Z();
   if(p1dotp2>0.) 
   {
      pass=false;
   }
   return pass;
}

event driver::GetEvent(double entry)
{
   GetEntry(entry);
   event evnt(*mu_pos_vec,*mu_neg_vec);
   return evnt;
}

double driver::GetEntries()
{
   return fChain->GetEntriesFast();
}


#endif // #ifdef driver_cxx
