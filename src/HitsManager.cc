#include "HitsManager.h"

const float TWID            =70.;
const int   N10MAX          =70;
const int   NhitsMAX        =90;
const int   N10THR          =6;

HitsManager::HitsManager() :
fTWindowLow( 6000. /* ns */ ),
fTWindowUp( 200000. /* ns */ ),
fRemoveDN( 0 ),
fRemoveAP( 0 )
{
    fAllHits=new HitCluster();
    this->ClearHits();
}

HitsManager::~HitsManager()
{
    if( fAllHits ){ delete fAllHits; fAllHits=NULL; }
    this->ClearHits();
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
    bool save=true;
    if( !(tim>fTWindowLow && tim<fTWindowUp) ){ save=false; }
    if( fRemoveDN && flg==kNoiseDN )          { save=false; }
    if( fRemoveAP && flg==kNoiseAP )          { save=false; }

    if( save )
    {
        fTimRaw.push_back(tim);
        fCabRaw.push_back(cab);
        fFlgRaw.push_back(flg);
    }
}

void HitsManager::SelectDigiHits()
{
    // Sort the selected PMT hits in order of incleasing time
    int nSelHits=int( fTimRaw.size() );
    std::vector<int> idx(nSelHits, -1);
    TMath::Sort(nSelHits, &fTimRaw[0], &idx[0], kFALSE);
    
    // Store the sorted PMT hits to the private memebers
    fAllHits->Clear();
    fAllHits->ReserveNhits(nSelHits);
    for(int i=0; i<nSelHits; i++)
    {
        fAllHits->AddHit(fCabRaw[idx[i]], fTimRaw[idx[i]], fFlgRaw[idx[i]]);
    }
    this->ClearHits();
}

void HitsManager::ClearHits()
{
    fTimRaw.clear();
    fCabRaw.clear();
    fFlgRaw.clear();
    fClusters.clear();
    std::vector<float>().swap(fTimRaw);
    std::vector<int>().swap(fCabRaw);
    std::vector<int>().swap(fFlgRaw);
    std::vector<HitCluster>().swap(fClusters);
}

void HitsManager::AddCluster(const int &first_hit, const int &last_hit, const int &n10)
{
    HitCluster aCluster;
    int flag=0;
    for(int j=first_hit; j<=last_hit; j++)
    {
        flag=0;
        if( fAllHits->S(j)!=0 ){ flag=-1; } // Non-Cherenkov hits
        aCluster.AddHit(fAllHits->C(j), 
                        fAllHits->T(j), 
                        flag);
    }
    aCluster.SetN10( n10 );
    fClusters.push_back( aCluster );
}

void HitsManager::SearchNeuCandidates()
{
    int     N10=0;
    int     nHitsCand=0;
    // !!! The initialization is important !!!
    float   prevT0=0.;
    int     prevIdx0=-999999;
    int     nCand=0;
    int     last_hit=-9999;

    this->SelectDigiHits();
    for(unsigned int i=0; i<fAllHits->N(); i++)
    {
        if( nCand>=1 && fAllHits->T(i)-fAllHits->T(prevIdx0)<TWID )
        {
            continue; 
        }

        N10=this->GetNXXX(fAllHits, i, 10., last_hit);
        if( N10>=N10THR )
        {
            if( prevIdx0<0 ) // The first candidate of this event
            {
                prevT0      =TWID + 1.;
            }
            else
            {
                prevT0      =fAllHits->T(i) - fAllHits->T(prevIdx0);
            }
            prevIdx0    =i;
            last_hit    =-9999;
            nHitsCand   =this->GetNXXX(fAllHits, i, 60., last_hit);
            if( prevT0>TWID && 
                N10<N10MAX  &&
                nHitsCand<NhitsMAX)
            {
                nCand+=1;
                this->AddCluster(i, last_hit, N10);
            }
        }
    }
}

int HitsManager::GetNXXX(const HitCluster *cl,
                         const int &start_idx, 
                         const float &twid,
                         int &last_idx)
{
    int n=0;
    float tdiff=0.;
    int nmax=cl->N();
    for(int i=start_idx; i<nmax; i++)
    {
        tdiff=cl->T(i) - cl->T(start_idx);
        if( i>nmax || tdiff>twid ){ last_idx=i-1; break; }
        else                      { n+=1; }
    }
    return n;
}
