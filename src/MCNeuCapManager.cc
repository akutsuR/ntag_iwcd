#include "MCNeuCapManager.h"

MCNeuCapManager::MCNeuCapManager()
{
}

MCNeuCapManager::~MCNeuCapManager()
{
        std::cout<<" Destructing" <<std::endl;
    std::vector<MCNeuCap_t>().swap(fNeuCap);
    std::cout<<" Hello: " <<std::endl;
}

void MCNeuCapManager::GetLineAgeInfo(const int &trgtIdx, 
                                     const vector<Track_t> &track, 
                                     int &grandID, 
                                     int &grandPDG, 
                                     bool &isPriNeu,
                                     int &ancestorPDG)
{
    int prntIdx=-999999;
    for(unsigned int i=0; i<track.size(); i++)
    {
        if( track[trgtIdx].parentID==track[i].thisID )
        {
            prntIdx=i; 
            break;
        }
    }

    grandID =track[prntIdx].thisID;
    grandPDG=track[prntIdx].thisPDG;
    isPriNeu=false; 
    if( grandPDG==GL::kPDG_NEUTRON &&
        track[prntIdx].parentID==0 )
    {
        isPriNeu=true;
        ancestorPDG=track[prntIdx].thisPDG;
    }
    else
    {
        int idxThis     =prntIdx;
        int ParentID    =-1;
        int idx         =-1;
        int nTry        =0;
        while( ParentID!=0 )
        {
            GetParentIndexAndTrackID(idxThis, track, ParentID, idx);
            if( idx<0 )
            {
                std::cout<<" [WARNING] MCNeuCapManager::GetLineAgeInfo "<<std::endl;
                std::cout<<"   - Cannot find parent particle for track ID: " << track[idxThis].thisID <<std::endl;
                break;
            }
            //std::cout<<" nTry: " << nTry 
            //         <<" - idxThis: " << idxThis
            //         <<" - ParetnID[idxThis]: " << track[idxThis].parentID
            //         <<" - idxNext: " << idx
            //         <<" - ParentID[idxNext]: " << track[idx].parentID
            //         <<" - ParentID: " << ParentID
            //         <<std::endl;
            idxThis=idx;
            nTry+=1;
            if( nTry>10 ){ exit(-1); }
        }
        if( idx>=0 ){ ancestorPDG=track[idx].thisPDG; }
        else        { ancestorPDG=-9999999; }
    }
}

void MCNeuCapManager::GetParentIndexAndTrackID(const int &trgtIdx, 
                                               const vector<Track_t> &track,
                                               int &prntId, 
                                               int &prntIdx)
{
    prntIdx=-999999;
    for(unsigned int i=0; i<track.size(); i++)
    {
        if( track[i].thisID==track[trgtIdx].parentID )
        {
            prntIdx=i; 
            prntId=track[i].parentID;
            break;
        }
    }
}


void MCNeuCapManager::FindMCTruthNCaptures(vector<Track_t> track)
{
    fNeuCap.clear();
    std::vector<MCNeuCap_t>().swap(fNeuCap);

    int index=-99999;
    int id=-99999;
    for(unsigned int i=0; i<track.size(); i++)
    {
        this->GetParentIndexAndTrackID(i, track, id, index);
        if( index<0 ){ continue; }
        if( track[i].thisPDG>GL::kPDG_NUCLEUS_OFFSET && 
            track[index].thisPDG==GL::kPDG_NEUTRON)
        {
            MCNeuCap_t aNCap;
            aNCap.nucleusPDG    =track[i].thisPDG-GL::kPDG_NUCLEUS_OFFSET;
            aNCap.captureTime   =track[i].time;
            aNCap.parentID      =track[i].parentID;
            aNCap.capturePos[0] =track[i].posSrt[0];
            aNCap.capturePos[1] =track[i].posSrt[1];
            aNCap.capturePos[2] =track[i].posSrt[2];

            int grandID     =-1;
            int grandPDG    =-1;
            bool isPriNeu   =false;
            int ancestorPDG =-1;
            GetLineAgeInfo(i, track, grandID, grandPDG, isPriNeu, ancestorPDG);
            aNCap.grandID       =grandID;
            aNCap.grandPDG      =grandPDG;
            aNCap.isPriNeu      =isPriNeu;
            aNCap.ancestorPDG   =ancestorPDG;
            fNeuCap.push_back( aNCap );
        }
    }
    if( fNeuCap.size()>0 )
    {
        this->FindGammas(track);
    }
}

void MCNeuCapManager::FindGammas(const vector<Track_t> &track)
{
    std::vector<bool> count(track.size(), false);
    for(unsigned int i=0; i<fNeuCap.size(); i++)
    {
        fNeuCap[i].nGammas=0;
        fNeuCap[i].eGammaTotal=0.;
        for(unsigned int j=0; j<track.size(); j++)
        {
            if( track[j].parentID==fNeuCap[i].parentID &&
                track[j].thisPDG==GL::kPDG_GAMMA &&
                count[j]==false )
            {
                fNeuCap[i].eGamma[fNeuCap[i].nGammas]=track[j].energy;
                fNeuCap[i].eGammaTotal+=track[j].energy;
                fNeuCap[i].nGammas+=1;
                count[j]=true; // Avoid double-counting
            }
        }
    }
}

double MCNeuCapManager::GetCaptureTime(const int &incap) const
{
    if( incap>=int( fNeuCap.size() ) )
    {
        std::cout<<" Current event has no " << incap <<" th capture" <<std::endl;
        std::cout<<" -> EXIT" <<std::endl;
        exit(-1);
    }
    return fNeuCap[incap].captureTime;
}

int MCNeuCapManager::GetCaptureNucleus(const int &incap) const
{
    if( incap>=int( fNeuCap.size() ) )
    {
        std::cout<<" Current event has no " << incap <<" th capture" <<std::endl;
        std::cout<<" -> EXIT" <<std::endl;
        exit(-1);
    }
    return fNeuCap[incap].nucleusPDG;
}

int MCNeuCapManager::GetNumGammas(const int &incap) const
{
    if( incap>=int( fNeuCap.size() ) )
    {
        std::cout<<" Current event has no " << incap <<" th capture" <<std::endl;
        std::cout<<" -> EXIT" <<std::endl;
        exit(-1);
    }
    return fNeuCap[incap].nGammas;
}

double MCNeuCapManager::GetTotalGammaEnergy(const int &incap) const
{
    if( incap>=int( fNeuCap.size() ) )
    {
        std::cout<<" Current event has no " << incap <<" th capture" <<std::endl;
        std::cout<<" -> EXIT" <<std::endl;
        exit(-1);
    }
    return fNeuCap[incap].eGammaTotal;
}




double MCNeuCapManager::GetCapturePosition(const int &incap, const int &ipos) const
{
    if( incap>=int( fNeuCap.size() ) )
    {
        std::cout<<" Current event has no " << incap <<" th capture" <<std::endl;
        std::cout<<" -> EXIT" <<std::endl;
        exit(-1);
    }
    return fNeuCap[incap].capturePos[ipos];
}
