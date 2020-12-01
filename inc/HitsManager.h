#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "TTree.h"
#include "TMath.h"
#include "TFile.h"

#include "global.h"
#include "HitCluster.h"

const int MAXHITS=200000;
const int NOISE_DN=1;
const int NOISE_AP=2;

class HitsManager
{
    public:
        HitsManager();
        virtual ~HitsManager();

        void ClearHits();
        void ReserveDigiHits(const int &);
        void AddDigiHit(const int&, const float&, const int&);
        void SelectDigiHits();
        void SearchNeuCandidates();

        int  GetNumOfNCandidates() const;
        int GetN10(const int&) const;
        float GetFirstHitTime(const int&) const;
        std::vector<HitCluster> GetClusters() const;

        // Only for MC studies
        void SetRemoveAfterpulse(const bool);
        void SetRemoveDarkNoise(const bool);

    private:
        void ReserveSelHits(const int&);
        void AddCluster(const int&);
        int GetNXXX(const int&, const float&);

        TTree *fDigiHits;
        Int_t fNDigiHits;
        Float_t fT[MAXHITS];
        Int_t fTube[MAXHITS];
        Bool_t fIsAfterpulse[MAXHITS];
        Bool_t fIsDark[MAXHITS];

        std::vector<float>  fTimRaw;
        std::vector<int>    fCabRaw;
        std::vector<int>    fFlgRaw;

        std::vector<float>  fTim;
        std::vector<int>    fCab;
        std::vector<int>    fFlg;
        std::vector<HitCluster> fClusters;

        float fTWindowLow;
        float fTWindowUp; 
        bool  fRemoveDN;
        bool  fRemoveAP;

        //////////////////////////
        Int_t fNCand;
        Float_t fCapTime[100];
        Int_t   fNXXX[7][100];
};
