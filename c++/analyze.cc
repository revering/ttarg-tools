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
   std::string fname, scalefile, strmap;
   double map;
   if(argc==4)
   {
      fname = argv[1];
      map = atof(argv[2]);
      scalefile = argv[3];
      strmap = argv[2];
   }
   else
   {
      std::cerr << "Please enter a root file, A' mass, and scaling file.\n";
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
   driver dver(itree, map, scalefile);
   std::string ofname = "Zll_dbrem_map_" + strmap + ".root";
   TFile *ofile = new TFile(ofname.c_str(),"RECREATE");
   TDirectory* inv_mass = ofile->mkdir("inv_mass_cut");
   TDirectory* after_dbrem = inv_mass->mkdir("after_dbrem");
   TDirectory* after_cuts = after_dbrem->mkdir("after_cuts");

   std::string all_name = "all_events";
   std::string inv_name = "inv_mass_cut";
   std::string cut_name = "after_cuts";
   std::string brem_name = "after_brem";
   dver.all_events.book(ofile,all_name);
   dver.after_cuts.book(after_cuts,cut_name);
   dver.after_brem.book(after_dbrem,brem_name);
   dver.inv_mass_cut.book(inv_mass,inv_name);
   TLorentzVector brem;
   event evnt, aft_brem;  
   int whichbrem = 0;
   int nprecuts = 0;
   int npostmass = 0;
   int npostpt = 0;
   int nposteta = 0;

   //histogram the initial events here
   double nentries = dver.GetEntries();
   for (Long64_t entry=0; entry<nentries; entry++)
   {
      evnt = dver.GetEvent(entry);
      dver.all_events.fill(evnt);
      nprecuts++;
      if(dver.pass_inv_mass(evnt)==true)
      {
         dver.inv_mass_cut.fill(evnt);
         npostmass++;
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
         if(whichbrem==1)
	 {
	    if((abs(evnt.mu_pos_vec.Eta())<2.4)&&(abs(evnt.mu_neg_vec.Eta())<2.1))
	    {
	       nposteta++;
	    }
	 }
	 if(whichbrem==-1)   
	 {
            if((abs(evnt.mu_pos_vec.Eta())<2.1)&&(abs(evnt.mu_neg_vec.Eta())<2.4))
	    {
	       nposteta++;
	    }
	 }
         //make cuts
         if (dver.pass_post_cuts(aft_brem)==true)
         {
            dver.after_cuts.fill(aft_brem);
	    npostpt++;
         }
      }	 
   }

   printf("Pre cuts: %d, Post mass: %d, Post eta: %d, Post pt: %d\n", nprecuts, npostmass, nposteta, npostpt);

//   dver.all_events.write();
//   dver.after_cuts.write(after_cuts);
 //  dver.after_brem.write(after_dbrem);
   
   ofile->Write();
   ofile->Close();
   f->Close();

   return 0;
}
