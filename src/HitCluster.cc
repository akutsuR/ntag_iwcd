#include "HitCluster.h"

HitCluster::HitCluster()
{
    this->Clear();
}

HitCluster::~HitCluster()
{
    this->Clear();
}

int HitCluster::N() const
{
    return int( fTime.size() );
}

int HitCluster::S(const int iHit) const
{
    return fSig[iHit];
}


void HitCluster::AddHit(const int &pmt_id,
                        const float &hit_time,
                        const int &sig_hit_flag)
{
    fPMTid.push_back( pmt_id );
    fTime.push_back( hit_time );
    fSig.push_back( sig_hit_flag );
}

void HitCluster::SetPreFitVertex(const float *vtx)
{
    fPreVtx[0]=vtx[0];
    fPreVtx[1]=vtx[1];
    fPreVtx[2]=vtx[2];
    fPreVtx[3]=vtx[3];
}

float HitCluster::GetTrueHitFraction() const
{
    int nSig=0;
    int nTot=int( fSig.size() );
    for(int i=0; i<nTot; i++)
    {
        if( fSig[i]==0 ){ nSig+=1; }
    }
    std::cout<<std::endl;
    return (float)nSig/nTot;
}


void HitCluster::Print()
{
    std::cout<<" ----- Pre-fit vertex: (" << fPreVtx[0]
             <<", " << fPreVtx[1]
             <<", " << fPreVtx[2]
             <<", " << fPreVtx[3]
             <<") "
             <<std::endl;
}
//------ Private ---------//
void HitCluster::Clear()
{
    fPMTid.clear();
    fTime.clear();
    fSig.clear();
    fPMTid.shrink_to_fit();
    fTime.shrink_to_fit();
    fSig.shrink_to_fit();

    fPreVtx[0]=-9999.;
    fPreVtx[1]=-9999.;
    fPreVtx[2]=-9999.;
    fPreVtx[3]=-9999.;
}

