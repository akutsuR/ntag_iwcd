// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "ResTFcn.h"
#include "VertexFit.h"

#include <cassert>

namespace ROOT {
   namespace Minuit2 {

ResTFcn::ResTFcn(VertexFit* vf) : 
fVF( vf ),
fDoTToFSDFit( false ),
fDoTToFT0Fit( false ),
fDoLikeliFit( false ),
fErrorDef(1.) 
{
}

double ResTFcn::operator()(const std::vector<double>& par) const
{
    assert( par.size()==4 || par.size()==3 );
    
    double vtx[4]={0.};
    if( par.size()==4 ){ for(int i=0; i<4; i++){ vtx[i]=par[i]; } }
    if( par.size()==3 ){ for(int i=0; i<3; i++){ vtx[i]=par[i]; } }
    
    for(int i=0; i<int( par.size() ); i++)
    {
        vtx[i]=par[i];
    }

    double metric=0.;
    if( fDoTToFSDFit ){ metric=fVF->GetTToFSD( vtx ); }
    if( fDoTToFT0Fit ){ metric=fVF->GetTToFT0( vtx ); }
    if( fDoLikeliFit ){ metric=fVF->GetNLLike( vtx ); }
    return metric;
}

void ResTFcn::PrintCurrentFlags()
{
    std::cout<<" DoTToFSDFit:" << fDoTToFSDFit <<std::endl;
    std::cout<<" DoTToFT0Fit:" << fDoTToFT0Fit <<std::endl;
    std::cout<<" DoLikeliFit:" << fDoLikeliFit <<std::endl;
}

  }  // namespace Minuit2
}  // namespace ROOT
