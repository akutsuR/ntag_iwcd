// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef MN_ResTFcn_H_
#define MN_ResTFcn_H_

#include <iostream>
#include "Minuit2/FCNBase.h"

class VertexFit;

namespace ROOT{
   namespace Minuit2{

class ResTFcn : public FCNBase
{
    public:
        ResTFcn(VertexFit* vf);
        virtual double Up() const {return fErrorDef;}
        virtual double operator()(const std::vector<double>&) const;

        void SetFlagLikeliFit(const bool doFit){ fDoLikeliFit=doFit; }    
        void SetFlagTToFT0Fit(const bool doFit){ fDoTToFT0Fit=doFit; }    
        void SetFlagTToFSDFit(const bool doFit){ fDoTToFSDFit=doFit; }    

        void PrintCurrentFlags();
    private:
        VertexFit *fVF;
        double fErrorDef;
        bool fDoTToFSDFit; // Fit by minimizing SD of t-tof
        bool fDoTToFT0Fit; // Fit by minimizing Gaussian t-tof-**t0**
        bool fDoLikeliFit; // Fit by maximizing a likelihood as a function of t-t0-tof
};

  }  // namespace Minuit2
}  // namespace ROOT
#endif //MN_GaussFcn_H_
