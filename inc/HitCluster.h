#pragma once
#include <iostream>
#include <string>
#include <vector>

class HitCluster
{
    public :
        HitCluster();
        virtual ~HitCluster();

        int N() const {return (int)fTim.size(); }
        int S(const int &i) const { return fFlg[i]; }
        int C(const int &i) const { return fCab[i]; }
        float T(const int &i) const { return fTim[i]; }
        int N10() const {return fN10; }

        void ReserveNhits(const int&);
        void AddHit(const int&, const float&, const int&);
        void SetN10(const int& n){ fN10=n; }
        void Clear();

        float GetTrueHitFraction() const;

    private :
        std::vector<int> fCab;
        std::vector<float> fTim;
        std::vector<int> fFlg;
        int fN10;
};
