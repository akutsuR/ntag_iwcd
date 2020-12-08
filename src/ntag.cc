#include <iostream>
#include <string>
#include <vector>
#include <time.h>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"

#include "WCSimRootGeom.hh"
#include "WCSimRootEvent.hh"
#include "TNRooTrackerVtx.hh"

#include "NeuCapMCInfo.h"
#include "MCNeuCapManager.h"
#include "HitsManager.h"
#include "VertexFit.h"
#include "NtagData.h"
#include "NtagUtil.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;


TTree *fRooTrackerTree=NULL;
TClonesArray *fnRooTrackerVtxTCA=NULL;
Double_t fDetCentreY=0.;


void SetThisNeuCandInfo(NtagData *da, const HitCluster &cl, const float *vtx)
{
    NtagUtil *nu=NtagUtil::GetInstance(); 
    int i=da->fncrNCand;
    da->fncrN10Raw[i]   =cl.N10();
    for(int j=0; j<4; j++)
    {
        da->fncrVtx[i][j]=vtx[j];
    }
    da->fncrDwall[i]    =nu->GetDwall(vtx);
    da->fncrTRaw0[i]    =cl.T(0);
    da->fncrNhits[i]    =cl.N();
    for(int j=0; j<cl.N(); j++)
    {
        da->fncrCabHits[i][j]   =cl.C(j);
        da->fncrTimHits[i][j]   =cl.T(j);
        da->fncrFlgHits[i][j]   =cl.S(j);
    }
    da->fncrNCand+=1;
}

void SetMCTruthMichelEInfo(NtagData *da, const MCNeuCapManager *nm)
{
    da->fnctNMiE=0;
    for(int i=0; i<nm->GetNumOfMichelEs(); i++)
    {
        for(int j=0; j<4; j++)
        {
            da->fnctMiEVtx[i][j]=nm->GetMiEVertex(i,j);
        }
        da->fnctMiEDwall[i]     =nm->GetMiEDwall(i);
        da->fnctMiEEnergy[i]    =nm->GetMiEEnergy(i);
        da->fnctMiEPrntPDG[i]   =nm->GetMiEParentPDG(i);
        da->fnctMiEGrPrntID[i]  =nm->GetMiEGrandParentID(i);
        da->fnctNMiE            +=1;
    }
}

void SetMCTruthNeuCapInfo(NtagData *da, const MCNeuCapManager *nm)
{
    da->fnctNCap=0;
    for(int i=0; i<nm->GetNumOfNeuCaptures(); i++)
    {
        for(int j=0; j<4; j++)
        {
            da->fnctVtx[i][j]=nm->GetCaptureVertex(i,j);
        }
        da->fnctDwall[i]=nm->GetCaptureDwall(i);
        da->fnctNuc[i]  =nm->GetCaptureNucleus(i);
        da->fnctType[i] =nm->GetAncestorType(i);
        da->fnctNGam[i] =nm->GetNumGammas(i);
        da->fnctETot[i] =nm->GetTotalGammaEnergy(i);
        da->fnctNCap    +=1;

        da->fnctNAntr[i]=nm->GetNumOfAncestors(i);
        for(int j=0; j<da->fnctNAntr[i]; j++)
        {
            da->fnctAntrEnergy[i][j]    =nm->GetAntrEnergy(i,j);
            da->fnctAntrVtx_x[i][j]     =nm->GetAntrVtx(i,j,0);
            da->fnctAntrVtx_y[i][j]     =nm->GetAntrVtx(i,j,1);
            da->fnctAntrVtx_z[i][j]     =nm->GetAntrVtx(i,j,2);
            da->fnctAntrVtx_t[i][j]     =nm->GetAntrVtx(i,j,3);
            da->fnctAntrDir_x[i][j]     =nm->GetAntrDir(i,j,0);
            da->fnctAntrDir_y[i][j]     =nm->GetAntrDir(i,j,1);
            da->fnctAntrDir_z[i][j]     =nm->GetAntrDir(i,j,2);
            da->fnctAntrPrntID[i][j]    =nm->GetAntrParentID(i,j);
            da->fnctAntrPDG[i][j]       =nm->GetAntrPDG(i,j);
        }
    }
}

void SetDigiHits(const WCSimRootTrigger *ev, HitsManager *hitMan, int &trigType)
{
    float toffset=950.; // Tentative
    trigType=0;
    if( ev->GetTriggerType()!=0 ){ toffset=100.; trigType=1;}

    int NDigiHits0=ev->GetNcherenkovdigihits();
    hitMan->ReserveDigiHits(NDigiHits0);
    for(int jHit=0; jHit<NDigiHits0; jHit++)
    {
        WCSimRootCherenkovDigiHit *aDigiHit= 
            dynamic_cast<WCSimRootCherenkovDigiHit*>(ev->GetCherenkovDigiHits()->At(jHit));
        float rawt=aDigiHit->GetT() - toffset;
        int   cabid=aDigiHit->GetTubeId();
        int   flag =0;  // Not implemented yet
        hitMan->AddDigiHit(cabid, rawt, 0);
    }
}

void GetTracks(const WCSimRootTrigger *ev, vector<Track_t> &trackArr)
{
    trackArr.clear();
    std::vector<Track_t>().swap(trackArr);
    trackArr.reserve(ev->GetNtrack()-2);
    for(int i=2; i<ev->GetNtrack(); i++)
    {
        WCSimRootTrack *wcsimroottrack 
                = dynamic_cast<WCSimRootTrack*>(ev->GetTracks()->At(i) );
        Track_t aTrack;
        aTrack.thisPDG      =wcsimroottrack->GetIpnu();
        aTrack.thisID       =wcsimroottrack->GetId();
        aTrack.parentID     =wcsimroottrack->GetParentId();
        aTrack.time         =wcsimroottrack->GetTime();
        aTrack.posSrt[0]    =wcsimroottrack->GetStart(0);
        aTrack.posSrt[1]    =wcsimroottrack->GetStart(1) - fDetCentreY;
        aTrack.posSrt[2]    =wcsimroottrack->GetStart(2);
        aTrack.energy       =wcsimroottrack->GetE();
        trackArr.push_back( aTrack );
    }
}

void PrintMCTruthNeuCap(const MCNeuCapManager *ncapMan)
{
    if( ncapMan->GetNumOfNeuCaptures()>0 )
    {
        for(int j=0; j<ncapMan->GetNumOfNeuCaptures(); j++)
        {
            cout<<" j: " << j
                <<" - time: " << ncapMan->GetCaptureVertex(j, 3)
                <<" - nuc: "  << ncapMan->GetCaptureNucleus(j)
                <<" - nG: " << ncapMan->GetNumGammas(j)
                <<" - eGTot: " << ncapMan->GetTotalGammaEnergy(j)
                <<" - pos: (" << ncapMan->GetCaptureVertex(j, 0)
                <<", " << ncapMan->GetCaptureVertex(j, 1)
                <<", " << ncapMan->GetCaptureVertex(j, 2)
                <<") - dwall: " << ncapMan->GetCaptureDwall(j)
                //<<" - gPDG: " << ncapMan->GetGrandParentPDG(j)
                //<<" - aPDG: " << ncapMan->GetAncestorParentPDG(j)
                //<<" - isPriNeu: " << ncapMan->GetIsPrimaryNeu(j)
                <<endl;
                if( j+1==ncapMan->GetNumOfNeuCaptures() ){ cout<<endl; }
        }
    }
}

void GetGeometryInfo(const TString filename)
{
    TFile *f=TFile::Open(filename);
    TTree *t=(TTree*)f->Get("wcsimGeoT");
    WCSimRootGeom *wcgeom=new WCSimRootGeom();
    t->SetBranchAddress("wcsimrootgeom", &wcgeom);
    t->GetEntry( 0 );

    fDetCentreY=wcgeom->GetWCOffset(1);
    cout<<" DetCentreY: " << fDetCentreY <<endl;
    NtagUtil *nu=NtagUtil::GetInstance();
    cout<<" Number of PMTs: " << wcgeom->GetWCNumPMT() <<endl;
    for(int i=0;  i<wcgeom->GetWCNumPMT(); i++)
    {
        WCSimRootPMT aPMT=wcgeom->GetPMT(i);
        TVector3 pos;
        TVector3 dir;
        pos(0)=aPMT.GetPosition(0);
        pos(1)=aPMT.GetPosition(1) - fDetCentreY;
        pos(2)=aPMT.GetPosition(2);
        dir(0)=aPMT.GetOrientation(0);
        dir(1)=aPMT.GetOrientation(1);
        dir(2)=aPMT.GetOrientation(2),
        nu->AddPMTInfo(pos, dir);
    }
    if( wcgeom ){ delete wcgeom;    wcgeom=NULL; }
    cout<<endl;
}

//ooooooooO000000000000000Ooooooooooo//
int main(int argc, char **argv)
{
    TString InFileName  =argv[1];
    TString OutFileName =argv[2];
///////////////////////////////////////////////

    cout<<" InFileName: " << InFileName <<endl;
    cout<<" OutFileName: " << OutFileName <<endl;

    GetGeometryInfo(InFileName);

    VertexFit *vf=new VertexFit();
    NtagUtil *nu=NtagUtil::GetInstance();
    vf->SetPMTPositions(nu->GetPMTPositions());

    TFile *file = TFile::Open(InFileName);
    //gSystem->Load("$WCSIMDIR/libWCSimRoot.so");
    //
    // Get the a pointer to the tree from the file
    TTree *tree = (TTree*)file->Get("wcsimT");
    fRooTrackerTree=(TTree*)file->Get("fRooTrackerOutputTree");
    fnRooTrackerVtxTCA=new TClonesArray("NRooTrackerVtx");
    fRooTrackerTree->SetBranchAddress("NRooTrackerVtx", &fnRooTrackerVtxTCA);
 
    // Create a WCSimRootEvent to put stuff from the tree
    WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();

    // Set the branch address for reading from the tree
    TBranch *branch = tree->GetBranch("wcsimrootevent");
    branch->SetAddress(&wcsimrootsuperevent);
    tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

    int nEntries=tree->GetEntries();
    vector<Track_t> trackArr;
    MCNeuCapManager *ncapMan=new MCNeuCapManager();
    HitsManager *hitMan=new HitsManager();

    NtagData *nData=new NtagData();
    nData->CreateTree(OutFileName);

    int trigType=-1;
    int nCandTot=0;
    int nNeuCapTot=0;
    clock_t start = clock();
    for(int iEntry=0; iEntry<nEntries; iEntry++)
    {
        tree->GetEntry(iEntry);      
        fRooTrackerTree->GetEntry( iEntry );
        NRooTrackerVtx *aVtx=(NRooTrackerVtx*)fnRooTrackerVtxTCA->At(0);

        WCSimRootTrigger *wcsimrootevent=NULL;
        wcsimrootevent=wcsimrootsuperevent->GetTrigger(0);
        
        GetTracks(wcsimrootevent, trackArr);
        ncapMan->FindMCTruthNCaptures( trackArr );
        ncapMan->FindMCTruthMichelEs( trackArr );
        //PrintMCTruthNeuCap(ncapMan);
        SetMCTruthNeuCapInfo(nData, ncapMan);
        SetMCTruthMichelEInfo(nData, ncapMan);

        SetDigiHits(wcsimrootevent, hitMan, trigType); 
        hitMan->SearchNeuCandidates();
        cout<<" iEntry: " << iEntry 
            <<" - mode: " << aVtx->EvtCode->GetString().Atoi()
            <<" - nMCNeu: " << ncapMan->GetNumOfNeuCaptures()
            <<" - nNeuCand: " << hitMan->GetNumOfNCandidates() 
            <<endl;
        cout<<endl;

        nCandTot+=hitMan->GetNumOfNCandidates();
        nNeuCapTot+=ncapMan->GetNumOfNeuCaptures();

        float vtxTmp[4]={0.};
        for(int k=0; k<hitMan->GetNumOfNCandidates(); k++)
        {
            HitCluster cl=hitMan->GetCluster(k);
            vf->FitVertex(&cl);
            for(int l=0; l<4; l++){ vtxTmp[l]=vf->GetNLLikeVertex(l); }
            SetThisNeuCandInfo(nData, cl, vtxTmp);
        }
        nData->fTrigType=trigType;
        nData->FillTree();
        wcsimrootsuperevent->ReInitialize();
        aVtx=NULL;
        fnRooTrackerVtxTCA->Clear();
    }
    clock_t end = clock();
    std::cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
    cout<<" nCandTot: " <<  nCandTot <<endl;

    nData->WriteTree();

return 0;
}
