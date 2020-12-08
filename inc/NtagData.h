#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "global.h"

using std::string;
using std::cout;
using std::endl;
using namespace GL;


class NtagData
{
    public:
        NtagData();
        virtual ~NtagData(){};
        void CreateTree(const TString&);
        void WriteTree();
        void SetTree(TTree*);
        void FillTree();

    private:
        void ClearVariables();
        void TryToMathCandWithTruth(); 
        void MatchToTrueNeuCap();
        void MatchToTrueMichelE();

        TFile *fFile;
        TTree *fTree;

    public:
        Int_t       fTrigType;

        // Reco. n capture including false event 
        Int_t       fncrNCand;
        Int_t       fncrN10Raw[kMaxRecoNCap];
        Float_t     fncrVtx[kMaxRecoNCap][4];
        Float_t     fncrDwall[kMaxRecoNCap];
        Float_t     fncrTRaw0[kMaxRecoNCap]; 
        Int_t       fncrIdxNCap[kMaxRecoNCap];
        Float_t     fncrDvtxNCap[kMaxRecoNCap][4];
        Int_t       fncrIdxMiE[kMaxRecoNCap];
        Float_t     fncrDtMiE[kMaxRecoNCap];
        Int_t       fncrNhits[kMaxRecoNCap];
        Int_t       fncrCabHits[kMaxRecoNCap][kNhitsMAX];
        Float_t     fncrTimHits[kMaxRecoNCap][kNhitsMAX];
        Int_t       fncrFlgHits[kMaxRecoNCap][kNhitsMAX];

        // MC truth n capture
        Int_t       fnctNCap;
        Float_t     fnctVtx[kMaxTrueNCap][4];
        Float_t     fnctDwall[kMaxTrueNCap];
        Int_t       fnctNuc[kMaxTrueNCap];
        Int_t       fnctType[kMaxTrueNCap];
        Int_t       fnctNGam[kMaxTrueNCap];
        Float_t     fnctETot[kMaxTrueNCap];

        // Ancestors associated with neutron capture
        Int_t       fnctNAntr[kMaxTrueNCap]; 
        Float_t     fnctAntrVtx_x[kMaxTrueNCap][kMaxTrueAntr];
        Float_t     fnctAntrVtx_y[kMaxTrueNCap][kMaxTrueAntr];
        Float_t     fnctAntrVtx_z[kMaxTrueNCap][kMaxTrueAntr];
        Float_t     fnctAntrVtx_t[kMaxTrueNCap][kMaxTrueAntr];
        Float_t     fnctAntrDir_x[kMaxTrueNCap][kMaxTrueAntr];
        Float_t     fnctAntrDir_y[kMaxTrueNCap][kMaxTrueAntr];
        Float_t     fnctAntrDir_z[kMaxTrueNCap][kMaxTrueAntr];
        Float_t     fnctAntrEnergy[kMaxTrueNCap][kMaxTrueAntr];
        Int_t       fnctAntrPrntID[kMaxTrueNCap][kMaxTrueAntr];
        Int_t       fnctAntrPDG[kMaxTrueNCap][kMaxTrueAntr];


        // MC truth Michel-e (to avoid confusion of n capture)
        Int_t       fnctNMiE;
        Float_t     fnctMiEVtx[kMaxTrueMiE][4];
        Float_t     fnctMiEDwall[kMaxTrueMiE];
        Float_t     fnctMiEEnergy[kMaxTrueMiE];
        Int_t       fnctMiEPrntPDG[kMaxTrueMiE];
        Int_t       fnctMiEGrPrntID[kMaxTrueMiE];
};
