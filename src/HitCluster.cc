#include "HitCluster.h"

HitCluster::HitCluster()
{
    this->Clear();
}

HitCluster::~HitCluster()
{
    this->Clear();
}


void HitCluster::AddHit(const int &pmt_id,
                        const float &hit_time,
                        const int &sig_hit_flag)
{
    fCab.push_back( pmt_id );
    fTim.push_back( hit_time );
    fFlg.push_back( sig_hit_flag );
}

float HitCluster::GetTrueHitFraction() const
{
    int nSig=0;
    int nTot=int( fFlg.size() );
    for(int i=0; i<nTot; i++)
    {
        if( fFlg[i]==0 ){ nSig+=1; }
    }
    std::cout<<std::endl;
    return (float)nSig/nTot;
}

void HitCluster::Clear()
{
    fCab.clear();
    fTim.clear();
    fFlg.clear();
    std::vector<float>().swap(fTim);
    std::vector<int>().swap(fFlg);
    std::vector<int>().swap(fCab);
}

void HitCluster::ReserveNhits(const int &n)
{
    fCab.reserve(n);
    fTim.reserve(n);
    fFlg.reserve(n);
}
