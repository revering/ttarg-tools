//To compile : g++ dbrem.cc -o simulate `root-config --cflags --glibs` DarkPhotons.cc -lgsl -lgslcblas -lm
#include "DarkPhotons.hh"
#include "driver.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
   std::string fname;
   double map;
   if(argc==3)
   {
      fname = argv[1];
      map = atof(argv[2]);
   }
   else
   {
      std::cerr << "Please enter a file name.\n";
      exit(1);
   }
   TFile *f = new TFile(fname.c_str(),"READ");
   if (!f || !f->IsOpen())
   {
      std::cerr << "Unable to open input root file.\n";
      exit(1);
   }
   TDirectory * indirec = f->GetDirectory("DirectoryName");
   TTree * itree = (TTree*)f->Get("all_events");
   driver dver(itree, map);
   TFile *ofile = new TFile("outputtest.root","Create");  
   TDirectory* after_cuts = f->mkdir("after_cuts");
   after_cuts->cd();
   TDirectory* after_dbrem = f->mkdir("after_dbrem");
   after_dbrem->cd();

   TLorentzVector brem;
   event evnt, aft_brem;  
   int whichbrem = 0;
   //histogram the initial events here
   double nentries = dver.GetEntries();
   for (Long64_t entry=0; entry<nentries; entry++)
   {
      evnt = dver.GetEvent(entry);
      dver.all_events.fill(evnt);

   //make cuts
      if (dver.cut(evnt)==false)
      {
         dver.after_cuts.fill(evnt);

         //do dark brem
         whichbrem = dver.Select_Lepton(evnt); 
         if (whichbrem == 1)
	 {
            brem = dver.simulate_dbrem(evnt.mu_pos_vec);
            aft_brem = event(brem, evnt.mu_neg_vec, whichbrem, evnt.mu_pos_vec);
         }
	 if (whichbrem == -1)
	 {
            brem = dver.simulate_dbrem(evnt.mu_neg_vec);
            aft_brem = event(evnt.mu_pos_vec, brem, whichbrem, evnt.mu_neg_vec);
         }
         dver.after_brem.fill(aft_brem);
       }   
   }

   //dver.all_events.write(itree);
   dver.after_cuts.write(after_cuts);
   dver.after_brem.write(after_dbrem);

   f->Close();

   return 0;
}
