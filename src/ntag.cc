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

void SetDigiHits(const WCSimRootTrigger *ev, HitsManager *hitMan)
{
    float toffset=950.; // Tentative
    if( ev->GetTriggerType()!=0 ){ toffset=100.; }

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
                <<" - time: " << ncapMan->GetCaptureTime(j)
                <<" - nuc: "  << ncapMan->GetCaptureNucleus(j)
                <<" - nG: " << ncapMan->GetNumGammas(j)
                <<" - eGTot: " << ncapMan->GetTotalGammaEnergy(j)
                <<" - pos: (" << ncapMan->GetCapturePosition(j, 0)
                <<", " << ncapMan->GetCapturePosition(j, 1)
                <<", " << ncapMan->GetCapturePosition(j, 2)
                <<") - gPDG: " << ncapMan->GetGrandParentPDG(j)
                <<" - aPDG: " << ncapMan->GetAncestorParentPDG(j)
                <<" - isPriNeu: " << ncapMan->GetIsPrimaryNeu(j)
                <<endl;
                if( j+1==ncapMan->GetNumOfNeuCaptures() ){ cout<<endl; }
        }
    }
}


//void GetGeometryInfo(const TString filename, VertexFit* vf)
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
}

//ooooooooO000000000000000Ooooooooooo//
int main(int argc, char **argv)
{
    TString InFileName  =argv[1];
    TString OutFileName =argv[2];
///////////////////////////////////////////////

    cout<<" InFileName: " << InFileName <<endl;
    cout<<" OutFileName: " << OutFileName <<endl;

    VertexFit *vf=new VertexFit();

    //GetGeometryInfo(InFileName, vf);
    GetGeometryInfo(InFileName);
    NtagUtil *nu=NtagUtil::GetInstance();
    vf->SetPMTPositions(nu->GetPMTPositions());

    TFile *file = TFile::Open(InFileName);
    gSystem->Load("$WCSIMDIR/libWCSimRoot.so");
    // Get the a pointer to the tree from the file
    TTree *tree = (TTree*)file->Get("wcsimT");
    fRooTrackerTree=(TTree*)file->Get("fRooTrackerOutputTree");
    fnRooTrackerVtxTCA=new TClonesArray("NRooTrackerVtx");
    fRooTrackerTree->SetBranchAddress("NRooTrackerVtx", &fnRooTrackerVtxTCA);
 
    // Create a WCSimRootEvent to put stuff from the tree in
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

    int nCandTot=0;
    int nNeuCapTot=0;
    clock_t start = clock();
    //for(int iEntry=0; iEntry<nEntries; iEntry++)
    for(int iEntry=0; iEntry<100; iEntry++)
    {
        tree->GetEntry(iEntry);      
        fRooTrackerTree->GetEntry( iEntry );
        NRooTrackerVtx *aVtx=(NRooTrackerVtx*)fnRooTrackerVtxTCA->At(0);

        WCSimRootTrigger *wcsimrootevent=NULL;
        wcsimrootevent=wcsimrootsuperevent->GetTrigger(0);
        
        GetTracks(wcsimrootevent, trackArr);
        ncapMan->FindMCTruthNCaptures( trackArr );
        //PrintMCTruthNeuCap(ncapMan);
        SetDigiHits(wcsimrootevent, hitMan); 
        hitMan->SearchNeuCandidates();

        cout<<" iEntry: " << iEntry 
            <<" - nMCNeu: " << ncapMan->GetNumOfNeuCaptures()
            <<" - nNeuCand: " << hitMan->GetNumOfNCandidates() 
            <<endl;
        cout<<endl;

        nCandTot+=hitMan->GetNumOfNCandidates();
        nNeuCapTot+=ncapMan->GetNumOfNeuCaptures();

        vector<HitCluster> cl=hitMan->GetClusters();
        for(int k=0; k<hitMan->GetNumOfNCandidates(); k++)
        {
            vf->SetCluster(&cl[k]);
            vf->FitVertex();
            //cout<<"    - " << k 
            //    <<" th cand - N10: " << hitMan->GetN10(k)
            //    <<" - tFirst: "      << hitMan->GetFirstHitTime(k)
            //    <<" - tRec: "        << vf->GetNLLikeVertex(3)
            //    <<" - Pos: ("        << vf->GetNLLikeVertex(0)
            //    <<", "        << vf->GetNLLikeVertex(1)
            //    <<", "        << vf->GetNLLikeVertex(2)
            //    <<")"
            //    <<endl;
        }
        //cout<<endl;
        //
        nData->FillTree();
        wcsimrootsuperevent->ReInitialize();
        aVtx=NULL;
        fnRooTrackerVtxTCA->Clear();
    }

    clock_t end = clock();
    std::cout << "duration = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
    cout<<" nCandTot: " <<  nCandTot <<endl;

return 0;
}
