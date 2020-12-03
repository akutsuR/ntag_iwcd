#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "global.h"
#include "NeuCapMCInfo.h"
#include "NtagUtil.h"

using std::vector;

class MCNeuCapManager
{
    public :
        MCNeuCapManager();
        virtual ~MCNeuCapManager();

        void FindMCTruthNCaptures(const vector<Track_t>&);
        void FindMCTruthMichelEs(const vector<Track_t>&);

        int GetNumOfNeuCaptures() const                             { return int(fNeuCap.size()); }
        int GetCaptureNucleus(const int &i) const                   { return fNeuCap[i].nucleusPDG; };
        float GetCaptureVertex(const int &i, const int &j) const    { return j<=2 ? fNeuCap[i].capturePos[j] : fNeuCap[i].captureTime; }
        float GetCaptureDwall(const int &i) const                   { return fNeuCap[i].captureDwall; }
        int GetNumGammas(const int &i) const                        { return fNeuCap[i].nGammas; }
        float GetTotalGammaEnergy(const int &i) const               { return fNeuCap[i].eGammaTotal; }
        int GetGrandParentPDG(const int &i) const                   { return fNeuCap[i].grandPDG; }
        int GetAncestorParentPDG(const int &i) const                { return fNeuCap[i].ancestorPDG; }
        int GetIsPrimaryNeu(const int &i) const                     { return fNeuCap[i].isPriNeu; }

        int GetNumOfMichelEs() const                                { return int(fMichelE.size()); }
        float GetMiEVertex(const int &i, const int &j) const        { return j<=2 ? fMichelE[i].Pos[j] : fMichelE[i].Time; }
        float GetMiEDwall(const int &i) const                       { return fMichelE[i].Dwall; }
        float GetMiEEnergy(const int &i) const                      { return fMichelE[i].energy; }
        int GetMiEGrandParentID(const int &i) const                 { return fMichelE[i].grandID; }
        int GetMiEParentPDG(const int &i) const                     { return fMichelE[i].parentPDG; }
                 
    private :
        void GetParentIndexAndTrackID(const int&, const vector<Track_t>&, int&, int&);
        void GetLineAgeInfo(const int&, const vector<Track_t>&, int&, int&, bool &, int&);
        void FindGammas(const vector<Track_t>&);
        std::vector<MCNeuCap_t> fNeuCap;
        std::vector<MichelE_t> fMichelE;
};
