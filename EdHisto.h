   
#ifndef __EdHisto_h
#define __EdHisto_h

#include "TH1F.h"
#include "TRandom2.h"
#include "TMath.h"

class EdHisto : public TH1F {
 public:
 EdHisto() : TH1F() { ps_Random = 0;};
 EdHisto(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup) :  TH1F(name,title,nbinsx,xlow,xup) {} ;
 EdHisto(const char *name,const char *title,Int_t nbinsx,const Float_t  *xbins) :  TH1F(name,title,nbinsx,xbins) {} ;
 EdHisto(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins) :   TH1F(name,title,nbinsx,xbins) {} ;
 EdHisto(const TVectorF &v) :   TH1F(v) {};
 EdHisto(const TH1F &h1f) :  TH1F(h1f) {};
  // virtual ~EdHisto();
private:  
   TRandom2 *ps_Random;

public:
   Double_t GetEdRandom()   const;
   void SetRandom(TRandom2 *fRandom) ;
   TRandom2  *GetRandom2() {return ps_Random;};
 
};

#endif//__EdHisto_h

