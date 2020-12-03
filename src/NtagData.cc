#include "NtagData.h"

NtagData::NtagData()
{
    this->ClearVariables();
}

void NtagData::FillTree()
{
    this->TryToMathCandWithTruth();

    fFile->cd();
    fTree->Fill();
    this->ClearVariables();
}

void NtagData::ClearVariables()
{
    fTrigType=-1;

    fnctNCap=0;
    for(int i=0; i<kMaxTrueNCap; i++)
    {
        fnctDwall[i]    =0.;
        fnctNuc[i]      =0;
        fnctType[i]     =0;
        fnctNGam[i]     =0;
        fnctETot[i]     =0.;
        for(int j=0; j<4; j++)
        {
            fnctVtx[i][j]=0.;
        }
    }

    fncrNCand=0;
    for(int i=0; i<kMaxRecoNCap; i++)
    {
        fncrN10Raw[i]   =0;
        fncrDwall[i]    =0.;
        fncrNhits[i]    =0;
        fncrTRaw0[i]    =0.; 
        fncrIdxNCap[i]  =0;
        for(int j=0; j<4; j++)
        {
            fncrVtx[i][j]       =0.;
            fncrDvtxNCap[i][j]  =0.;
        }
        fncrIdxMiE[i]   =0;
        fncrDtMiE[i]    =0.; 
    }
}

void NtagData::CreateTree(const TString &OutFileName)
{
    fFile=new TFile(OutFileName, "recreate");
    fTree=new TTree("ntag", "");

    fTree->Branch("TrigType",       &fTrigType,     "TrigType/I");

    fTree->Branch("nctNCap",        &fnctNCap,          "nctNCap/I");
    fTree->Branch("nctVtx",         fnctVtx,            "nctVtx[nctNCap][4]/F");
    fTree->Branch("nctDwall",       fnctDwall,          "nctDwall[nctNCap]/F");
    fTree->Branch("nctNuc",         fnctNuc,            "nctNuc[nctNCap]/I");
    fTree->Branch("nctType",        fnctType,           "nctType[nctNCap]/I");
    fTree->Branch("nctNGam",        fnctNGam,           "nctNGam[nctNCap]/I");
    fTree->Branch("nctETot",        fnctETot,           "nctETot[nctNCap]/F");

    fTree->Branch("nctNMiE",        &fnctNMiE,          "nctNMiE/I");
    fTree->Branch("nctMiEVtx",      fnctMiEVtx,         "nctMiEVtx[nctNMiE][4]/F");
    fTree->Branch("nctMiEDwall",    fnctMiEDwall,       "nctMiEDwall[nctNMiE]/F");
    fTree->Branch("nctMiEEnergy",   fnctMiEEnergy,      "nctMiEEnergy[nctNMiE]/F");
    fTree->Branch("nctMiEPrntPDG",  fnctMiEPrntPDG,     "nctMiEPrntPDG[nctNMiE]/I");
    fTree->Branch("nctMiEGrPrntID", fnctMiEGrPrntID,    "nctMiEGrPrntID[nctNMiE]/I");

    fTree->Branch("ncrNCand",       &fncrNCand,      "ncrNCand/I");
    fTree->Branch("ncrN10Raw",      fncrN10Raw,      "ncrN10Raw[ncrNCand]/I");
    fTree->Branch("ncrVtx",         fncrVtx,         "ncrVtx[ncrNCand][4]/F");
    fTree->Branch("ncrDwall",       fncrDwall,       "ncrDwall[ncrNCand]/F");
    fTree->Branch("ncrNhits",       fncrNhits,       "ncrNhits[ncrNCand]/I");
    fTree->Branch("ncrTRaw0",       fncrTRaw0,       "ncrTRaw0[ncrNCand]/F");
    fTree->Branch("ncrIdxNCap",     fncrIdxNCap,     "ncrIdxNCap[ncrNCand]/I");
    fTree->Branch("ncrDvtxNCap",    fncrDvtxNCap,    "ncrDvtxNCap[ncrNCand][4]/F");
    fTree->Branch("ncrIdxMiE",      fncrIdxMiE,      "ncrIdxMiE[ncrNCand]/I");
    fTree->Branch("ncrDtMiE",       fncrDtMiE,       "ncrDtMiE[ncrNCand]/F");
}

void NtagData::SetTree(TTree *t)
{
    fTree=t;

    fTree->SetBranchAddress("TrigType",       &fTrigType);

    fTree->SetBranchAddress("nctNCap",        &fnctNCap);
    fTree->SetBranchAddress("nctVtx",         fnctVtx);
    fTree->SetBranchAddress("nctDwall",       fnctDwall);
    fTree->SetBranchAddress("nctNuc",         fnctNuc);
    fTree->SetBranchAddress("nctType",        fnctType);
    fTree->SetBranchAddress("nctNGam",        fnctNGam);
    fTree->SetBranchAddress("nctETot",        fnctETot);

    fTree->SetBranchAddress("nctNMiE",        &fnctNMiE);
    fTree->SetBranchAddress("nctMiEVtx",      fnctMiEVtx);
    fTree->SetBranchAddress("nctMiEDwall",    fnctMiEDwall);
    fTree->SetBranchAddress("nctMiEEnergy",   fnctMiEEnergy);
    fTree->SetBranchAddress("nctMiEPrntPDG",  fnctMiEPrntPDG);
    fTree->SetBranchAddress("nctMiEGrPrntID", fnctMiEGrPrntID);

    fTree->SetBranchAddress("ncrNCand",       &fncrNCand);
    fTree->SetBranchAddress("ncrN10Raw",      fncrN10Raw);
    fTree->SetBranchAddress("ncrVtx",         fncrVtx);
    fTree->SetBranchAddress("ncrDwall",       fncrDwall);
    fTree->SetBranchAddress("ncrNhits",       fncrNhits);
    fTree->SetBranchAddress("ncrTRaw0",       fncrTRaw0);
    fTree->SetBranchAddress("ncrIdxNCap",     fncrIdxNCap);
    fTree->SetBranchAddress("ncrDvtxNCap",    fncrDvtxNCap);
    fTree->SetBranchAddress("ncrIdxMiE",      fncrIdxMiE);
    fTree->SetBranchAddress("ncrDvtxMiE",     fncrDtMiE);
}

void NtagData::WriteTree()
{
    fFile->cd();
        fTree->Write();
    fFile->Close();
}



void NtagData::TryToMathCandWithTruth()
{
// - Actuall, do not perfrom matching each neuron candidate 
//   with true neutron capture
// - Find out the true neutron capture/Michel-e giving the shortest time difference 
//   from it, for each neutron capture/Michel-e instead.
// - By using this information, one can define whether or not 
//   each candidate corresponds to a MC truth neutron capture later 
//   (e.g. identify a candidate as a MC truth capture 
//   if its shortest time difference is <100ns)
// - Similary, one can check confusion of a candidate with Michel-e 
//   (The current definition of Michel-e includes all secondary e+/e-)

    this->MatchToTrueNeuCap();
    this->MatchToTrueMichelE();
}

void NtagData::MatchToTrueNeuCap()
{
    float aDt;
    float aDtMin=9999999.;
    int   idxMin=-99999;

    for(int i=0; i<fncrNCand; i++)
    {
        fncrIdxNCap[i]=-9999; 
        for(int j=0; j<4; j++){ fncrDvtxNCap[i][j]=-99999.; }
        
        aDtMin=999999999.;
        idxMin=-999999.;
        if( fnctNCap==0 ){ continue; }

        for(int j=0; j<fnctNCap; j++)
        {
            aDt=TMath::Abs(fnctVtx[j][3] - fncrVtx[i][3]);
            if( j==0 )
            { 
               aDtMin=aDt; 
               idxMin=j;
            }
            else
            {
                if( aDt<aDtMin )
                { 
                   aDtMin=aDt; 
                   idxMin=j;
                }
            }
        }
        if( aDtMin>399.9 ){ aDtMin=399.9; }

        for(int j=0; j<3; j++)
        {
            fncrDvtxNCap[i][j]=fncrVtx[i][j] - fnctVtx[idxMin][j]; 
        }
        fncrDvtxNCap[i][3]   =aDtMin; 
        fncrIdxNCap[i]       =idxMin;
    }
}

void NtagData::MatchToTrueMichelE()
{
    float aDt;
    float aDtMin=9999999.;
    int   idxMin=-99999;

    for(int i=0; i<fncrNCand; i++)
    {
        fncrIdxMiE[i]=-9999; 
        fncrDtMiE[i]=-999.;
        
        aDtMin=999999999.;
        idxMin=-999999.;
        if( fnctNMiE==0 ){ continue; }

        for(int j=0; j<fnctNCap; j++)
        {
            aDt=TMath::Abs(fnctMiEVtx[j][3] - fncrVtx[i][3]);
            if( j==0 )
            { 
               aDtMin=aDt; 
               idxMin=j;
            }
            else
            {
                if( aDt<aDtMin )
                { 
                   aDtMin=aDt; 
                   idxMin=j;
                }
            }
        }
        if( aDtMin>399.9 ){ aDtMin=399.9; }

        fncrDtMiE[i]      =aDtMin; 
        fncrIdxMiE[i]     =idxMin;
    }
}
