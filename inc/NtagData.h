#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"

using std::string;
using std::cout;
using std::endl;

const int kMaxTrueNCap=500;
const int kMaxRecoNCap=500;

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
        TFile *fFile;
        TTree *fTree;

    public:
        Int_t       fnctNCap;
        Float_t     fnctVtx[kMaxTrueNCap][4];
        Float_t     fnctDwall[kMaxTrueNCap];
        Int_t       fnctNuc[kMaxTrueNCap];
        Int_t       fnctType[kMaxTrueNCap];
        Int_t       fnctNGam[kMaxTrueNCap];
        Float_t     fnctETot[kMaxTrueNCap];

        Int_t       fncrNCand;
        Int_t       fncrN10Raw[kMaxRecoNCap];
        Float_t     fncrVtx[kMaxRecoNCap][4];
        Float_t     fncrDwall[kMaxRecoNCap];
        Int_t       fncrNhits[kMaxRecoNCap];
        Float_t     fncrTRaw0[kMaxRecoNCap]; 
        Int_t       fncrIdxMCT[kMaxRecoNCap];
        Float_t     fncrDvtxMCT[kMaxRecoNCap][4];
};
