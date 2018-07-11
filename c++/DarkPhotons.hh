#ifndef DarkPhotons_h
#define DarkPhotons_h 1

#define Mel 5.1E-4 // electron mass (GeV)
#define Mmu 0.1056 // muon mass (GeV)
#define alphaEW 1.0/137.0
#define MUp 2.79 // protonMu
#define Mpr 0.938 // proton mass
#define max_uint 4294967296.0l
#define GeVtoPb 3.894E+08

#define NPTAB 15
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "TLorentzVector.h"

struct ParamsForChi {double AA; double ZZ; double MMA; double EE0;};
struct momentum {double E0; double Theta; double Phi;};

class DarkPhotons
{

  public:

    DarkPhotons(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                double epsilIn, std::string fname);

    ~DarkPhotons();

    void SetSigmaNorm(double SigmaNormIn);
    double TotalCrossSectionCalc(double E0);
    double TotalMuCrossSectionCalc(double E0);
    double MaxCrossSectionCalc(double E0);
    void PrepareTable();
    double GetMA() {return MA;}
    double GetEThresh() {return EThresh;}
    double GetSigmaNorm() {return SigmaNorm;}
    double Getepsil() {return epsil;}
     // usage of normalization below:   Nsign = (Naccepted/Nsimulated)*Normalization*EOT
    double GetNormalization() {return 3.0e-15 * (Density/11.35) * (207./ANucl) *
                                      epsil * epsil / (SigmaNorm * epsilBench * epsilBench);}
    double GetsigmaTot(double E0);
    double GetsigmaMax(double E0);
    bool Emission(double E0, double DensityMat, double StepLength); // E0 in GeV, density in g/cm3, StepLength in mm
    momentum SimulateEmission(double E0);
    TLorentzVector* MuSimulateEmission(double E0);
    double GetAccumulatedProbability() {return AccumulatedProbability;}
    double GetMadgraphE(double E0);

  private:

    double MA;
    double EThresh;
    double SigmaNorm;
    double ANucl;
    double ZNucl;
    double Density;
    double epsilBench;
    double epsil;
    std::map< double , std::vector < double > > efracs;
    std::vector < std::pair < double, int > > energies;
    std::ifstream ifile;
    int nptable;
    double ep[15];
    double sigmap[15];
    double sigmax[15];

    double AccumulatedProbability;

};

#endif
