#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "global.h"
#include "NeuCapMCInfo.h"

using std::vector;

class MCNeuCapManager
{
    public :
        MCNeuCapManager();
        virtual ~MCNeuCapManager();

        //void FindMCTruthNCaptures(const vector<Track_t>&);
        void FindMCTruthNCaptures(vector<Track_t>);

        int GetNumOfNeuCaptures() const { return int(fNeuCap.size()); }
        double GetCaptureTime(const int&) const;
        int GetCaptureNucleus(const int&) const;
        double GetCapturePosition(const int&, const int&) const;
        int GetNumGammas(const int&) const;
        double GetTotalGammaEnergy(const int&) const;
        int GetGrandParentPDG(const int &j) const {return fNeuCap[j].grandPDG; }
        int GetAncestorParentPDG(const int &j) const {return fNeuCap[j].ancestorPDG; }
        int GetIsPrimaryNeu(const int &j) const {return fNeuCap[j].isPriNeu; }
    

    private :
        void GetParentIndexAndTrackID(const int&, const vector<Track_t>&, int&, int&);
        void GetLineAgeInfo(const int&, const vector<Track_t>&, int&, int&, bool &, int&);
        void FindGammas(const vector<Track_t>&);
        std::vector<MCNeuCap_t> fNeuCap;
};
