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
            idxThis=idx;
            nTry+=1;
            if( nTry>50 ){ exit(-1); }
        }
        if( idx>=0 ){ ancestorPDG=track[idx].thisPDG; }
        else        { ancestorPDG=-9999999; }
    }
}

bool MCNeuCapManager::GetParentIndexAndTrackID(const int &trgtIdx, 
                                               const vector<Track_t> &track,
                                               int &prntId, 
                                               int &prntIdx)
{
    prntIdx=-999999;
    bool isFound=false;
    for(unsigned int i=0; i<track.size(); i++)
    {
        if( track[i].thisID==track[trgtIdx].parentID )
        {
            prntIdx=i; 
            prntId=track[i].parentID;
            isFound=true;
            break;
        }
    }
    return isFound;
}


void MCNeuCapManager::FindMCTruthNCaptures(const vector<Track_t>& track)
{
    fNeuCap.clear();
    std::vector<MCNeuCap_t>().swap(fNeuCap);

    int prntIdx=-99999; // index of parent in ``track'' for given child
    int dummyID=-99999; // parent track ID of the parent
    NtagUtil *us=NtagUtil::GetInstance();
    std::vector<float> ctime;
    std::vector<MCNeuCap_t> NeuCapTmp;

    for(unsigned int i=0; i<track.size(); i++)
    {
        this->GetParentIndexAndTrackID(i, track, dummyID, prntIdx);
        if( prntIdx<0 ){ continue; }

        if( track[i].thisPDG>GL::kPDG_NUCLEUS_OFFSET && 
            track[prntIdx].thisPDG==GL::kPDG_NEUTRON)
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

            //int prntIdx=-1; 
            int strIdx=prntIdx;
            int n=0;
            bool isFound=false;
            while( 1 )
            {
                isFound=this->GetParentIndexAndTrackID(strIdx, track, dummyID, prntIdx);
                GetLineAgeInfo(i, track, grandID, grandPDG, isPriNeu, ancestorPDG);
                aNCap.antrVtx_x[n]  =track[strIdx].posSrt[0];
                aNCap.antrVtx_y[n]  =track[strIdx].posSrt[1];
                aNCap.antrVtx_z[n]  =track[strIdx].posSrt[2];
                aNCap.antrVtx_t[n]  =track[strIdx].time;
                aNCap.antrDir_x[n]  =track[strIdx].dirSrt[0];
                aNCap.antrDir_y[n]  =track[strIdx].dirSrt[1];
                aNCap.antrDir_z[n]  =track[strIdx].dirSrt[2];
                aNCap.antrEnergy[n] =track[strIdx].energy;
                aNCap.antrPrntID[n] =track[strIdx].parentID;
                aNCap.antrPDG[n]    =track[strIdx].thisPDG;
                n+=1;
                if( !isFound ){ break; }
                strIdx=prntIdx;
            }
            aNCap.nAncestor=n;

            //aNCap.grandID       =grandID;
            //aNCap.grandPDG      =grandPDG;
            //aNCap.isPriNeu      =isPriNeu;
            //aNCap.ancestorPDG   =ancestorPDG;

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
            // Cherenkov threshold, assuming
            // a refractive index of 1.34 for
            // water
            if( track[i].energy<0.768 ) 
            {
                cout<<" Found a Michel-e, but below Cherenkov threshold " <<endl;
            }
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


int MCNeuCapManager::GetAncestorType(const int &i) const
{
    // idex for the earliest ancestor
    int idxAntr0=fNeuCap[i].nAncestor-1;
    int type=-1;
    if( idxAntr0==0 && 
        fNeuCap[i].antrPDG[ idxAntr0 ]==GL::kPDG_NEUTRON &&
        fNeuCap[i].antrPrntID[ idxAntr0 ]==0 )
    {
        // Neutron after FSI is captured directory
        type=GL::eDirectPriNeu;
    }
    else if( fNeuCap[i].antrPDG[ idxAntr0 ]==GL::kPDG_PROTON &&
             fNeuCap[i].antrPrntID[ idxAntr0 ]==0 )
    {
        // Neutron produced via proton SI is captured
        type=GL::eProtonSecNeu;
    }
    else if( fNeuCap[i].antrPDG[ idxAntr0 ]==GL::kPDG_PiPlus &&
             fNeuCap[i].antrPrntID[ idxAntr0 ]==0 )
    {
        // Neutron produced via neutron SI is captured
        type=GL::ePiPlusSecNeu;
    }
    else if( fNeuCap[i].antrPDG[ idxAntr0 ]==GL::kPDG_PiMinus &&
             fNeuCap[i].antrPrntID[ idxAntr0 ]==0 )
    {
        // Neutron produced via pi+ SI is captured
        type=GL::ePiMinusSecNeu;
    }
    else if( fNeuCap[i].antrPDG[ idxAntr0 ]==GL::kPDG_PiMinus &&
             fNeuCap[i].antrPrntID[ idxAntr0 ]==0 )
    {
        // Neutron produced via pi- SI is captured
        type=GL::ePiMinusSecNeu;
    }

    else if( fNeuCap[i].antrPDG[ idxAntr0 ]==GL::kPDG_MuMinus &&
             fNeuCap[i].antrPrntID[ idxAntr0 ]==0 )
    {
        // Neutron produced via pi- SI is captured
        type=GL::eMuMinusCapNeu;
    }
    else
    {
        type=GL::eOtherSecNeu;
    }
    return type;
}
