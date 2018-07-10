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
      virtual ~event();
};

event::event(TLorentzVector pvec, TLorentzVector mvec, bool posbrem, TLorentzVector pre_vec)
{
   mu_pos_vec = pvec;
   mu_neg_vec = mvec;
   pbrem = posbrem;
   pre_brem_vec =  pre_vec;
}

event::event(TLorentzVector pvec, TLorentzVector mvec)
{
   mu_pos_vec = pvec;
   mu_neg_vec = mvec;
}

event::~event()
{
}


