#include "NtagData.h"

NtagData::NtagData()
{
    this->ClearVariables();
}

void NtagData::FillTree()
{
    fFile->cd();
    fTree->Fill();
    this->ClearVariables();
}

void NtagData::ClearVariables()
{
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
        fncrIdxMCT[i]   =0;
        for(int j=0; j<4; j++)
        {
            fncrVtx[i][j]       =0.;
            fncrDvtxMCT[i][j]   =0.;
        }
    }
}

void NtagData::CreateTree(const TString &OutFileName)
{
    fFile=new TFile(OutFileName, "recreate");
    fTree=new TTree("ntag", "");

    fTree->Branch("nctNCap",        &fnctNCap,      "nctNCap/I");
    fTree->Branch("nctVtx",         fnctVtx,        "nctVtx[nctNCap][4]/F");
    fTree->Branch("nctDwall",       fnctDwall,      "nctDwall[nctNCap]/F");
    fTree->Branch("nctNuc",         fnctNuc,        "nctNuc[nctNCap]/I");
    fTree->Branch("nctType",        fnctType,       "nctType[nctNCap]/I");
    fTree->Branch("nctNGam",        fnctNGam,       "nctNGam[nctNCap]/I");
    fTree->Branch("nctETot",        fnctETot,       "nctETot[nctNCap]/F");

    fTree->Branch("ncrNCand",       &fncrNCand,      "ncrNCand/I");
    fTree->Branch("ncrN10Raw",      fncrN10Raw,      "ncrN10Raw[ncrNCand]/I");
    fTree->Branch("ncrVtx",         fncrVtx,         "ncrVtx[ncrNCand][4]/F");
    fTree->Branch("ncrDwall",       fncrDwall,       "ncrDwall[ncrNCand]/F");
    fTree->Branch("ncrNhits",       fncrNhits,       "ncrNhits[ncrNCand]/I");
    fTree->Branch("ncrTRaw0",       fncrTRaw0,       "ncrTRaw0[ncrNCand]/F");
    fTree->Branch("ncrIdxMCT",      fncrIdxMCT,      "ncrIdxMCT[ncrNCand]/I");
    fTree->Branch("ncrDvtxMCT",     fncrDvtxMCT,     "ncrDvtxMCT[ncrNCand][4]/F");
}

void NtagData::SetTree(TTree *t)
{
    fTree=t;
    fTree->SetBranchAddress("nctNCap",        &fnctNCap);
    fTree->SetBranchAddress("nctVtx",         fnctVtx);
    fTree->SetBranchAddress("nctDwall",       fnctDwall);
    fTree->SetBranchAddress("nctNuc",         fnctNuc);
    fTree->SetBranchAddress("nctType",        fnctType);
    fTree->SetBranchAddress("nctNGam",        fnctNGam);
    fTree->SetBranchAddress("nctETot",        fnctETot);

    fTree->SetBranchAddress("ncrNCand",       &fncrNCand);
    fTree->SetBranchAddress("ncrN10Raw",      fncrN10Raw);
    fTree->SetBranchAddress("ncrVtx",         fncrVtx);
    fTree->SetBranchAddress("ncrDwall",       fncrDwall);
    fTree->SetBranchAddress("ncrNhits",       fncrNhits);
    fTree->SetBranchAddress("ncrTRaw0",       fncrTRaw0);
    fTree->SetBranchAddress("ncrIdxMCT",      fncrIdxMCT);
    fTree->SetBranchAddress("ncrDvtxMCT",     fncrDvtxMCT);
}

void NtagData::WriteTree()
{
    fFile->cd();
        fTree->Write();
    fFile->Close();
}
