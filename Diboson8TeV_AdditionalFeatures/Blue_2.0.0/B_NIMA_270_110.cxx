//-----------------------------------------------------------------------------
//
// BLUE: A ROOT class implementing the Best Linear Unbiased Estimate method.
// 
// Copyright (C) 2012-2014, Richard.Nisius@mpp.mpg.de
// All rights reserved
//
// This file is part of BLUE - Version 2.0.0.
//
// BLUE is free software: you can redistribute it and/or modify it under the 
// terms of the GNU Lesser General Public License as published by the Free 
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// For the licensing terms see the file COPYING or http://www.gnu.org/licenses.
//
//-----------------------------------------------------------------------------
#include "Blue.h"

void B_NIMA_270_110(){
  //----------------------------------------------------------------------------
  // Only one result, Equation 19, is calculated
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst = 4;
  static const Int_t NumUnc = 1;
  static const Int_t NumObs = 1;

  // Index for which estimate determines which observable
  Int_t IWhichObs[NumEst] = {0, 0, 0, 0};

  // Input is based on the results Eq(14)
  // and the uncertainty matrix Eq(17')
  // 2.74, 1.15, 0.86, 1.31,
  // 1.15, 1.67, 0.82, 1.32,
  // 0.86, 0.82, 2.12, 1.05,
  // 1.31, 1.32, 1.05, 2.93
  // 1) get sigma_i = sqrt(cov_ii)
  // 2) get rho_ij from cov_ij=rho_ij/sigma_i/sigma_j

  // Fill estimates for single uncertainty
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEstOne[LenXEst] = {
     9.5, 1.66,
    11.9, 1.29,
    11.1, 1.45,
     8.9, 1.71
  };

  static const Int_t LenCor = NumEst * NumEst;
  Double_t Cor[LenCor] = {
    1.00, 0.54, 0.36, 0.46,
    0.54, 1.00, 0.44, 0.60,
    0.36, 0.44, 1.00, 0.42,
    0.46, 0.60, 0.42, 1.00
  };

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->PrintStatus();

  // Fill all estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEstOne[ind]);
    ind = ind + NumUnc + 1;
  }

  // Fill correlation
  for(Int_t k = 0; k<NumUnc; k++){
      myBlue->FillCor(k, &Cor[0]);
  }
  
  // Fix input and inspect input
  myBlue->FixInp();
  myBlue->PrintEst();
  myBlue->PrintCor(0);
  myBlue->PrintCov();
  myBlue->PrintRho();

  // Solve and inspect result
  myBlue->Solve();
  printf("... B_NIMA_270_110: Calculate Equation %2i \n", 19);
  myBlue->PrintResult();

  // Delete Object
  delete myBlue;
  myBlue = NULL;
  return;
}
