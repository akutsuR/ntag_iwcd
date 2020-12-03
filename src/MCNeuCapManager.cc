#include "MCNeuCapManager.h"

MCNeuCapManager::MCNeuCapManager()
{
}

MCNeuCapManager::~MCNeuCapManager()
{
    std::vector<MCNeuCap_t>().swap(fNeuCap);
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


void MCNeuCapManager::FindMCTruthNCaptures(const vector<Track_t>& track)
{
    fNeuCap.clear();
    std::vector<MCNeuCap_t>().swap(fNeuCap);

    int index=-99999;
    int id=-99999;
    NtagUtil *us=NtagUtil::GetInstance();
    std::vector<float> ctime;
    std::vector<MCNeuCap_t> NeuCapTmp;

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
            aNCap.captureDwall  =us->GetDwall(track[i].posSrt);

            int grandID     =-1;
            int grandPDG    =-1;
            bool isPriNeu   =false;
            int ancestorPDG =-1;
            GetLineAgeInfo(i, track, grandID, grandPDG, isPriNeu, ancestorPDG);
            aNCap.grandID       =grandID;
            aNCap.grandPDG      =grandPDG;
            aNCap.isPriNeu      =isPriNeu;
            aNCap.ancestorPDG   =ancestorPDG;
            NeuCapTmp.push_back( aNCap );

            ctime.push_back( track[i].time );
        }
    }

    int nNeuCap=(int)ctime.size();
    if( nNeuCap )
    {
        std::vector<int> idx(nNeuCap, -1);
        TMath::Sort(nNeuCap, &ctime[0], &idx[0], kFALSE);
        fNeuCap.reserve(nNeuCap);   
        for(int i=0; i<nNeuCap; i++){ fNeuCap.push_back( NeuCapTmp[ idx[i] ] ); }

        if( fNeuCap.size()>0 )
        {
            this->FindGammas(track);
        }
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

void MCNeuCapManager::FindMCTruthMichelEs(const vector<Track_t> &track)
{
    fMichelE.clear();
    std::vector<MichelE_t>().swap(fMichelE);

    NtagUtil *us=NtagUtil::GetInstance();
    int id=-99999;
    int index=-9999;
    for(unsigned int i=0; i<track.size(); i++)
    {
        if( track[i].parentID==0 ){ continue; }

        this->GetParentIndexAndTrackID(i, track, id, index);
        if( index<0 ){ continue; }
        if( track[i].thisPDG==GL::kPDG_ELECTRON )
        {
            cout<<" Found a Michel-e" <<endl;
            MichelE_t aMiE;
            aMiE.Time       =track[i].time;
            aMiE.Pos[0]     =track[i].posSrt[0];
            aMiE.Pos[1]     =track[i].posSrt[1];
            aMiE.Pos[2]     =track[i].posSrt[2];
            aMiE.Dwall      =us->GetDwall(track[i].posSrt);
            aMiE.energy     =track[i].energy;
            aMiE.grandID    =track[index].parentID;
            aMiE.parentPDG  =track[index].thisPDG; 
            fMichelE.push_back( aMiE );
        }
    }
}
