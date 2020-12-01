#pragma once
#include <iostream>
#include <string>
#include <vector>

class HitCluster
{
    public :
        HitCluster();
        virtual ~HitCluster();

        int N() const;  
        int C(const int &i) const { return fPMTid[i]; }
        float T(const int &i) const { return fTime[i]; }
        int S(const int) const;
        float GetTrueHitFraction() const;

        void AddHit(const int&, const float&, const int&);
        void SetPreFitVertex(const float*);
        void Print();

    private :
        void Clear();
        std::vector<int> fPMTid;
        std::vector<float> fTime;
        std::vector<int> fSig;
        float fPreVtx[4];
};
