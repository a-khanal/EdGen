#include "EdHisto.h"


Double_t EdHisto::GetEdRandom() const
{
   // return a random number distributed according the histogram bin contents.
   // This function checks if the bins integral exists. If not, the integral
   // is evaluated, normalized to one.
   // The integral is automatically recomputed if the number of entries
   // is not the same then when the integral was computed.
   // NB Only valid for 1-d histograms. Use GetRandom2 or 3 otherwise.
   // If the histogram has a bin with negative content a NaN is returned 

   if (fDimension > 1) {
      Error("GetEdRandom","Function only valid for 1-d histograms");
      return 0;
   }
   Int_t nbinsx = GetNbinsX();
   Double_t integral = 0;
   // compute integral checking that all bins have positive content (see ROOT-5894)
   if (fIntegral) {
      if (fIntegral[nbinsx+1] != fEntries) integral = ((TH1*)this)->ComputeIntegral(true);
      else  integral = fIntegral[nbinsx];
   } else {
      integral = ((TH1*)this)->ComputeIntegral(true);
   }
   if (integral == 0) return 0;
   // return a NaN in case some bins have negative content
   if (integral == TMath::QuietNaN() ) return TMath::QuietNaN(); 

   Double_t r1 = ps_Random->Rndm();
   Int_t ibin = TMath::BinarySearch(nbinsx,fIntegral,r1);
   Double_t x = GetBinLowEdge(ibin+1);
   if (r1 > fIntegral[ibin]) x +=
      GetBinWidth(ibin+1)*(r1-fIntegral[ibin])/(fIntegral[ibin+1] - fIntegral[ibin]);
   return x;
}


void EdHisto::SetRandom(TRandom2 *fRandom) {
  ps_Random = fRandom;
  printf("Set up Random Generator to Default TRandom2 seed=%d \n",ps_Random->GetSeed());

}
