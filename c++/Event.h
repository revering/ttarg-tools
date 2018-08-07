/**
 * Basic event class for storing Z-ll leptons.
 * Michael Revering, University of Minnesota.
 */

#ifndef event_h
#define event_h

#include "TLorentzVector.h"

class event
//Stores the positive and negative leptons, as well as the information of the lepton before 
//bremming. pbrem sets which lepton was selected, and is 1 for positive, -1 for negative, and 
//zero for neither.
{
   public:
      TLorentzVector mu_pos_vec;
      TLorentzVector mu_neg_vec;
      int           pbrem;
      TLorentzVector pre_brem_vec;

      event(TLorentzVector pvec, TLorentzVector mvec, int posbrem, TLorentzVector pre_vec);
      event(TLorentzVector pvec, TLorentzVector mvec);
      event();
      virtual ~event();
};

#endif
