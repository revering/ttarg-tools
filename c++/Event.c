#include "TLorentzVector.h"
#include "Event.h"

event::event()
{
}

event::event(TLorentzVector pvec, TLorentzVector mvec, bool posbrem, TLorentzVector pre_vec)
{
   mu_pos_vec = pvec;
   mu_neg_vec = mvec;
   pbrem = posbrem;
   pre_brem_vec = pre_vec;
}

event::event(TLorentzVector pvec, TLorentzVector mvec)
{
   mu_pos_vec = pvec;
   mu_neg_vec = mvec;
}

event::~event()
{
}


