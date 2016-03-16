// Based on TGenPhaseSpace
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Phase Space Generator, based on the GENBOD routine of CERNLIB           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
   
#ifndef __EdGenPhaseSpace_h
#define __EdGenPhaseSpace_h

#include "TLorentzVector.h"


class EdGenPhaseSpace : public TObject {
private:  
   Int_t        fNt;             // number of decay particles
   Double_t     fMass[18];       // masses of particles
   Double_t     fBeta[3];        // betas of decaying particle
   Double_t     fTeCmTm;         // total energy in the C.M. minus the total mass
   Double_t     fWtMax;          // maximum weigth 
   TLorentzVector  fDecPro[18];  //kinematics of the generated particles 

   Double_t PDK(Double_t a, Double_t b, Double_t c);  

public:
   EdGenPhaseSpace(): fNt(0), fMass(), fBeta(), fTeCmTm(0.), fWtMax(0.) {}
   EdGenPhaseSpace(const EdGenPhaseSpace &gen);
   virtual ~EdGenPhaseSpace() {}
   EdGenPhaseSpace& operator=(const EdGenPhaseSpace &gen);

   Bool_t          SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass, Option_t *opt="");
   Double_t        Generate();
   TLorentzVector *GetDecay(Int_t n); 

   Int_t    GetNt()      const { return fNt;}
   Double_t GetWtMax()   const { return fWtMax;}

};

#endif//__EdGenPhaseSpace_h

