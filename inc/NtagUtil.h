#pragma once

#include <iostream>
#include <vector>

#include "TVector3.h"
#include "TMath.h"

using std::cout;
using std::endl;
using std::vector;

class NtagUtil
{
    public :
        static NtagUtil* GetInstance();
        void Finalize();

        float GetDwall(const float*);
        float GetIDTankRadius(){ return fIDTankR; }
        float GetIDTankHeight(){ return fIDTankH; }

        vector<TVector3> GetPMTPositions() const {return fPos; }


        void AddPMTInfo(const TVector3&, const TVector3&);
        void SetIDTankRadius(const float &f){ fIDTankR=f; }
        void SetIDTankHalfHeight(const float &f){ fIDTankH=f; }

    protected :
        NtagUtil();

    private :
        static NtagUtil *fTheInstance;

        float fIDTankR;
        float fIDTankH;
        float fIDTankY;
        float fCWater;

        vector<TVector3> fPos; 
        vector<TVector3> fDir;
};
