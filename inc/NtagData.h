#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"

using std::string;
using std::cout;
using std::endl;

const int kMaxTrueNCap=500;
const int kMaxRecoNCap=500;
const int kMaxTrueMiE=100;

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

        // MC truth n capture
        Int_t       fnctNCap;
        Float_t     fnctVtx[kMaxTrueNCap][4];
        Float_t     fnctDwall[kMaxTrueNCap];
        Int_t       fnctNuc[kMaxTrueNCap];
        Int_t       fnctType[kMaxTrueNCap];
        Int_t       fnctNGam[kMaxTrueNCap];
        Float_t     fnctETot[kMaxTrueNCap];

        // MC truth Michel-e (to avoid confusion of n capture)
        Int_t       fnctNMiE;
        Float_t     fnctMiEVtx[kMaxTrueMiE][4];
        Float_t     fnctMiEDwall[kMaxTrueMiE];
        Float_t     fnctMiEEnergy[kMaxTrueMiE];
        Int_t       fnctMiEPrntPDG[kMaxTrueMiE];
        Int_t       fnctMiEGrPrntID[kMaxTrueMiE];

        // Reco. n capture including false event 
        Int_t       fncrNCand;
        Int_t       fncrN10Raw[kMaxRecoNCap];
        Float_t     fncrVtx[kMaxRecoNCap][4];
        Float_t     fncrDwall[kMaxRecoNCap];
        Int_t       fncrNhits[kMaxRecoNCap];
        Float_t     fncrTRaw0[kMaxRecoNCap]; 
        Int_t       fncrIdxNCap[kMaxRecoNCap];
        Float_t     fncrDvtxNCap[kMaxRecoNCap][4];
        Int_t       fncrIdxMiE[kMaxRecoNCap];
        Float_t     fncrDtMiE[kMaxRecoNCap];
};
