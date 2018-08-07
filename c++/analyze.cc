/**
 * @file analyze.cc
 * @brief Script performing dark brem simulation and analysis on Z->ll samples.
 * @author Michael Revering, University of Minnesota
 */

//To compile : g++ analyze.cc -o analyze `root-config --cflags --glibs` DarkPhotons.cc Event.c Plotter.c driver.c -lgsl -lgslcblas -lm
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
   //Read in the z->ll file, the A' mass, and the scaling file.
   std::string fname, scalefile, strmap;
   double map;
   if(argc==4)
   {
      fname = argv[1];
      map = atof(argv[2]);
      scalefile = argv[3];
      strmap = argv[2]; //Lazy way to get the A' mass as a string and an int.
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

   //Setup the directory structure for the output root file.
   TDirectory * indirec = f->GetDirectory("DirectoryName");
   TTree * itree = (TTree*)f->Get("all_events");
   driver dver(itree, map, scalefile);
   std::string ofname = "Zll_dbrem_map_" + strmap + ".root";
   TFile *ofile = new TFile(ofname.c_str(),"RECREATE");
   TDirectory* inv_mass = ofile->mkdir("inv_mass_cut");
   TDirectory* after_cuts = inv_mass->mkdir("after_cuts");
   TDirectory* after_dbrem = after_cuts->mkdir("after_dbrem");

   //Make the histograms in each of the subdirectories.
   std::string all_name = "all_events";
   std::string inv_name = "inv_mass_cut";
   std::string cut_name = "after_cuts";
   std::string brem_name = "after_brem";
   dver.all_events.book(ofile,all_name);
   dver.after_cuts.book(after_cuts,cut_name);
   dver.after_brem.book(after_dbrem,brem_name);
   dver.inv_mass_cut.book(inv_mass,inv_name);

   TLorentzVector brem;
   event evnt, aft_brem, aft_cuts;  
   int whichbrem = 0;
  
   //Variables to count numbers of events to pass each cut.
   int nprecuts = 0;
   int npostmass = 0;
   int npostpt = 0;
  
   //Loop over all events in the file.
   double nentries = dver.GetEntries();
   for (Long64_t entry=0; entry<nentries; entry++)
   {
      evnt = dver.GetEvent(entry);
      dver.all_events.fill(evnt); //Put each event into the "all events" histogram.
      nprecuts++;
      if(dver.pass_inv_mass(evnt)==true) //Make the invariant mass cut first.
      {
         dver.inv_mass_cut.fill(evnt);
         npostmass++;
	 //Choose which particle brems. Is weighted by the cross section at their respective energies. Returns 1 if the positive lepton brems, -1 if the negative brems.
         whichbrem = dver.Select_Lepton(evnt); 
         if (whichbrem == 1)
         {
	    //Do the brem. Returns a TLorentzVector of the particle after it brems.
            brem = dver.simulate_dbrem(evnt.mu_pos_vec); 
	    //Make a new event with the bremmed particle.
            aft_brem = event(brem, evnt.mu_neg_vec, whichbrem, evnt.mu_pos_vec);
	    //Make a new event with the old particle (needed so that I can plot the particle that bremmed before it changed, without changing any of the plotting code).
	    aft_cuts = event(evnt.mu_pos_vec, evnt.mu_neg_vec, whichbrem, evnt.mu_pos_vec);
         }
         else if (whichbrem == -1) //Same procedure for the negative one, if it was the brem.
         {
            brem = dver.simulate_dbrem(evnt.mu_neg_vec);
            aft_brem = event(evnt.mu_pos_vec, brem, whichbrem, evnt.mu_neg_vec);
            aft_cuts = event(evnt.mu_pos_vec, evnt.mu_neg_vec, whichbrem, evnt.mu_neg_vec);
	 }
         else {std::cout << "Skipping an event.\n";} //If both cross sections are zero.
	 
         if (dver.pass_cuts(aft_cuts)==true) //Make cuts now that we have a track and a muon.
         {
	    dver.after_cuts.fill(aft_cuts); //Plot histograms showing the track's information.
            dver.after_brem.fill(aft_brem); //Plot histograms showing how the track changed.
	    npostpt++;
         }
      }	 
   }

   printf("Pre cuts: %d, Post mass: %d, Post cuts: %d\n", nprecuts, npostmass, npostpt);

//   dver.all_events.write();
//"Write" method is currently used to make the integration histograms. The files are actually written in the last step, when we call ofile->Write().
   dver.after_cuts.write(); 
   dver.after_brem.write();
  
   ofile->Write();   
   ofile->Close();
   f->Close();

   return 0;
}
