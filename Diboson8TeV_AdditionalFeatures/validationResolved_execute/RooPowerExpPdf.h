/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOPOWEREXPPDF
#define ROOPOWEREXPPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooPowerExpPdf : public RooAbsPdf {
public:
  RooPowerExpPdf() {} ; 
  RooPowerExpPdf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _c,
	      RooAbsReal& _power);
  RooPowerExpPdf(const RooPowerExpPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooPowerExpPdf(*this,newname); }
  inline virtual ~RooPowerExpPdf() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy c ;
  RooRealProxy power ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooPowerExpPdf,1) // Your description goes here...
};
 
#endif