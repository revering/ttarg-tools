//To compile : g++ dbrem.cc -o simulate `root-config --cflags --glibs` DarkPhotons.cc -lgsl -lgslcblas -lm
#include "DarkPhotons.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
   double map, ebeam;
   if(argc==2)
   {
      std::string fname = argv[1];
   }
   std::string fname = argv[1];
   TFile *f = new TFile(fname.c_str(),"READ");
   TDirectory * indirec = f->GetDirectory("DirectoryName");
   TTree * itree = (TTree*)f->Get("treename");
   TLorentzVector * pvec = new TLorentzVector();
   TLorentzVector * mvec = new TLorentzVector();

   itree->SetBranchAddress("MuMinus",pvec);
   itree->SetBranchAddress("MuPlus",mvec);

   TDirectory* after_dbrem = f->mkdir("after_dbrem");
   after_dbrem->cd();
   TTree * otree = new TTree("Signal", "Tree containing Lorentz Vectors after Dark Brem");
   TLorentzVector * ovecp = new TLorentzVector();
   TLorentzVector * ovecm = new TLorentzVector();
   TLorentzVector * outvec = new TLorentzVector();
   bool * pbrem = new bool;
   otree->Branch("MuPlus","TLorentzVector",ovecp);
   otree->Branch("MuMinus","TLorentzVector",ovecm);
   otree->Branch("PreBrem","TLorentzVector",outvec);
   otree->Branch("WhichBrem","Bool",pbrem);

   std::string ifname = "/local/cms/user/revering/dphoton/code/geant/test.txt";
   DarkPhotons* dphoton = new DarkPhotons(map, 0, 1, 28, 14, 2.32, 1, ifname);
   dphoton->PrepareTable();
   double e_mom, a_z, a_t, e_z, e_t, e_x, e_y, a_E, a_y, a_x;
   int entries = itree->GetEntries();
   for(int i=0;i<entries;i++)
   {
      indirec->cd();
      itree->GetEntry(i);
      double pEnergy = pvec->E();
      double mEnergy = mvec->E();
      after_dbrem->cd();
      double pWeight = dphoton->GetsigmaTot(pEnergy);
      double mWeight = dphoton->GetsigmaTot(mEnergy);
      double A = drand48()*(pWeight+mWeight);
      if(A<pWeight)
      {
         outvec->SetPxPyPzE(pvec->X(),pvec->Y(),pvec->Z(),pvec->E());
         TVector3 UzVec = outvec->Vect().Unit();
         ovecp = (TLorentzVector*)dphoton->MuSimulateEmission(pEnergy)->Clone();
         ovecp->RotateUz(UzVec);
         ovecm = (TLorentzVector*)mvec->Clone();
	 *pbrem = true;
      }
      else
      {
         outvec->SetPxPyPzE(mvec->X(),mvec->Y(),mvec->Z(),mvec->E());
         TVector3 UzVec = outvec->Vect().Unit();
         ovecm = (TLorentzVector*)dphoton->MuSimulateEmission(mEnergy)->Clone();
         ovecm->RotateUz(UzVec);
         ovecp = (TLorentzVector*)pvec->Clone();
         *pbrem = false;
      }
     
      otree->Fill();
   }

   f->Write();
   f->Close();

   return 0;
}
