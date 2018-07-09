import ROOT

class Event:
   """
   Class for tracker as target events.
   """

   def __init__(self, positive_lepton, negative_lepton, pbrem = -1, after_brem = ROOT.TLorentzVector()):
      self.lepplus = lepton1
      self.lepminus = lepton2
      self.plusbrem = pbrem
      self.abrem = after_brem

