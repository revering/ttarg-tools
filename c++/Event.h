#ifndef event_h
#define event_h

#include "TLorentzVector.h"

class event
{
   public:
      TLorentzVector mu_pos_vec;
      TLorentzVector mu_neg_vec;
      bool           pbrem;
      TLorentzVector pre_brem_vec;

      event(TLorentzVector pvec, TLorentzVector mvec, bool posbrem, TLorentzVector pre_vec);
      event(TLorentzVector pvec, TLorentzVector mvec);
      event();
      virtual ~event();
};

#endif
