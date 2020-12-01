#include "NtagUtil.h"

NtagUtil* NtagUtil::fTheInstance=0;

NtagUtil::NtagUtil()
{
    fCWater =21.5833;   // m/ns
    fIDTankR=400.;
    fIDTankH=600.;
    
    fPos.clear();
    fDir.clear();
}

void NtagUtil::Finalize()
{
    vector<TVector3>().swap(fPos);
    vector<TVector3>().swap(fDir);
}

NtagUtil* NtagUtil::GetInstance()
{
    if( fTheInstance==0 )
    {
        fTheInstance=new NtagUtil();
        cout<<" Instanciating NtagUtil for the first time" <<endl;
    }
    //cout<<" Calling TheInstance " <<endl;
    return fTheInstance;
}

void NtagUtil::AddPMTInfo(const TVector3 &p, const TVector3 &d)
{
    fPos.push_back( p );
    fDir.push_back( d );
}
