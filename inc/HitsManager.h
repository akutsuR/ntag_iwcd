#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "TTree.h"
#include "TMath.h"
#include "TFile.h"

#include "global.h"
#include "HitCluster.h"

const int kNoiseDN=1;
const int kNoiseAP=2;

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

        int  GetNumOfNCandidates() const { return (int)fClusters.size(); }
        int GetN10(const int &i) const { return fClusters[i].N10(); }
        float GetFirstHitTime(const int &i) const { return fClusters[i].T(0); }
        HitCluster GetCluster(const int &i) const { return fClusters[i]; }

        std::vector<HitCluster> GetClusters() const { return fClusters; }

        // Only for MC studies
        void SetRemoveAfterpulse(const bool &b){ fRemoveAP=b; }
        void SetRemoveDarkNoise(const bool &b){ fRemoveDN=b; }

    private:
        void AddCluster(const int&, const int&, const int&);
        int GetNXXX(const HitCluster*, const int&, const float&, int&);

        std::vector<float>  fTimRaw;
        std::vector<int>    fCabRaw;
        std::vector<int>    fFlgRaw;
        std::vector<HitCluster> fClusters;
        HitCluster *fAllHits;

        float fTWindowLow;
        float fTWindowUp; 
        bool  fRemoveDN;
        bool  fRemoveAP;

        //////////////////////////
};
