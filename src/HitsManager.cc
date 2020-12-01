#include "HitsManager.h"

HitsManager::HitsManager() :
fTWindowLow( 6000. /* ns */ ),
fTWindowUp( 200000. /* ns */ ),
fRemoveDN( 0 ),
fRemoveAP( 0 )
{
    this->ClearHits();
}

HitsManager::~HitsManager()
{
    this->ClearHits();
}

void HitsManager::SetRemoveAfterpulse(const bool flag)
{
    fRemoveAP=flag;
}

void HitsManager::SetRemoveDarkNoise(const bool flag)
{
    fRemoveDN=flag;
}

void HitsManager::ReserveDigiHits(const int &n)
{
    fTimRaw.reserve(n);
    fCabRaw.reserve(n);
    fFlgRaw.reserve(n);
}


void HitsManager::AddDigiHit(const int &cab,
                             const float &tim,
                             const int &flg)
{
    fTimRaw.push_back(tim);
    fCabRaw.push_back(cab);
    fFlgRaw.push_back(flg);
}


void HitsManager::SelectDigiHits()
{
    const int nmax=int( fTimRaw.size() );
    std::vector<float> tim;
    std::vector<int> cab;
    std::vector<int> flg;
    tim.reserve(nmax);
    flg.reserve(nmax);
    cab.reserve(nmax);

    // Select PMT hits that are used for neutron tagging
    for(unsigned int i=0; i<fTimRaw.size(); i++)
    {
        if( fTimRaw[i]>fTWindowLow && fTimRaw[i]<fTWindowUp )
        {
            // Do not include dark noise
            if( fRemoveDN && fFlgRaw[i]==1 ){ continue; }

            // Do not include afterpulses
            if( fRemoveAP && fFlgRaw[i]==2 ){ continue; }

            tim.push_back( fTimRaw[i] );
            cab.push_back( fCabRaw[i] );
            flg.push_back( fFlgRaw[i] );
        }
    }

    // Sort the selected PMT hits in order of incleasing time
    int nSelHits=int( tim.size() );
    std::vector<int> idx(nSelHits, -1);
    TMath::Sort(nSelHits, &tim[0], &idx[0], kFALSE);
    
    // Store the sorted PMT hits to the private memebers
    this->ClearHits();
    this->ReserveSelHits(nSelHits);
    for(int i=0; i<nSelHits; i++)
    {
        fTim.push_back( tim[ idx[i] ] );
        fCab.push_back( cab[ idx[i] ] );
        fFlg.push_back( flg[ idx[i] ] );
    }
}


void HitsManager::ClearHits()
{
    fTimRaw.clear();
    fCabRaw.clear();
    fFlgRaw.clear();

    fTim.clear();
    fCab.clear();
    fFlg.clear();

    fClusters.clear();

    std::vector<float>().swap(fTimRaw);
    std::vector<int>().swap(fCabRaw);
    std::vector<int>().swap(fFlgRaw);
    std::vector<float>().swap(fTim);
    std::vector<int>().swap(fCab);
    std::vector<int>().swap(fFlg);
    std::vector<HitCluster>().swap(fClusters);
}

void HitsManager::ReserveSelHits(const int &nhits)
{
    fTim.reserve( nhits );
    fCab.reserve( nhits );
    fFlg.reserve( nhits );
}

void HitsManager::AddCluster(const int &first_hit)
{
    HitCluster aCluster;
    int j=0;
    int flag=0;
    for(int i=0; i<this->GetNXXX(first_hit, 60.); i++)
    {
        j=first_hit + i;
        flag=0;
        if( fFlg[j]!=0 ){ flag=-1; } // Non-Cherenkov hits
        aCluster.AddHit(fCab[j], fTim[j], flag);
    }
    fClusters.push_back( aCluster );
}

std::vector<HitCluster> HitsManager::GetClusters() const
{
    return fClusters;
}

void HitsManager::SearchNeuCandidates()
{
    const float TWID            =70.;
    const int   N10MAX          =70;
    const int   NhitsMAX        =90;
    const int   N10THR          =6;

    int     N10=0;
    int     nHitsCand=0;
    // !!! The initialization is important !!!
    float   prevT0=0.;
    int     prevIdx0=-999999;
    int     nCand=0;
    std::vector<float> ctim;

    const int N_NXXX_TYPES=7;
    const float TWID_NXXX[N_NXXX_TYPES]={
                    10.,
                    50.,
                    100.,
                    150., 
                    200.,
                    250.,
                    300.};
    std::vector< std::vector<int> > nxx;

    this->SelectDigiHits();
    for(int i=0; i<int( fTim.size() ); i++)
    {
        if( nCand>=1 && fTim[i]-fTim[prevIdx0]<TWID )
        {
            continue; 
        }

        N10=this->GetNXXX(i, 10.);
        if( N10>=N10THR )
        {
            if( prevIdx0<0 ) // The first candidate of this event
            {
                prevT0      =TWID + 1.;
            }
            else
            {
                prevT0      =fTim[i] - fTim[prevIdx0];
            }
            prevIdx0    =i;
            nHitsCand   =this->GetNXXX(i, 60.);
            if( prevT0>TWID && 
                N10<N10MAX  &&
                nHitsCand<NhitsMAX)
            {
                nCand+=1;
                ctim.push_back( fTim[i] );
                std::vector<int> inxx;
                for(int j=0; j<N_NXXX_TYPES; j++)
                {
                    inxx.push_back( this->GetNXXX(i, TWID_NXXX[j]) );
                }
                nxx.push_back( inxx );
                this->AddCluster(i);
            }
        }
    }

    fNCand=nCand;
    for(int i=0; i<int( ctim.size() ); i++)
    {
        fCapTime[i]=ctim[i];
    }

    for(int i=0; i<int( nxx.size() ); i++)
    {
        for(int j=0; j<N_NXXX_TYPES; j++)
        {
            fNXXX[j][i]=nxx[i][j];
        }
    }
}


float HitsManager::GetFirstHitTime(const int &icand) const
{
    return fCapTime[icand];
}


int HitsManager::GetN10(const int &icand) const
{
    return fNXXX[0][icand];
}


int HitsManager::GetNXXX(const int &start_idx, const float &twid)
{
    int n=0;
    float tdiff=0.;
    for(int i=start_idx; i<int( fTim.size() ); i++)
    {
        tdiff=fTim[i] - fTim[start_idx];
        if( i>int( fTim.size() ) || tdiff>twid )
        {
            break;
        }
        else
        {
            n+=1;
        }
    }
    return n;
}

int  HitsManager::GetNumOfNCandidates() const
{
    return fNCand;
}
