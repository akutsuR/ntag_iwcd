#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MinosError.h"


#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMath.h"
#include "TGraph.h"

#include "HitCluster.h"
#include "ResTFcn.h"
#include "NtagUtil.h"

using namespace ROOT::Minuit2;

const int kNAxis=4;
const int kNVtx=5;
const int kNEdge=2;

enum eAxis{ eX=0, eY, eZ, eT};
enum eVtx{ eTrue=0, eSeed=1, eTToFSD=2, eTToFT0=3, eNLLike=4};
enum eEdge{ eLow=0, eUp};
const vector<std::string> kSAxis{"x", "y", "z", "t"};

class VertexFit
{
    public:
        VertexFit();
        virtual ~VertexFit();

        void SetPMTPositions(const vector<TVector3> &p){ fPos=p; }
        void LoadTemplate(std::string);
        void SetCluster(HitCluster *aCl){ fCl=aCl; }
        void FitVertex();

        float GetNLLikeVertex(const int &i) const {return fVertex[eNLLike][i]; }
        float GetTToFT0Vertex(const int &i) const {return fVertex[eTToFT0][i]; }
        
        double GetTToFSD(const double*);
        double GetTToFT0(const double*);
        double GetNLLike(const double*);
        double GetLogLikeli(const double&);

        float GetTSD() const {return fTSD; }
        float GetAveTTOF() const {return fAveTTOF; }
        int GetMinPoint(const int &i) const{ return fMinPoint[i]; }

        void SetTOFFSET(const double &tofst){ fTOFFSET=tofst; }
        void SetTrueVertex(const double*);        


    private:
        TVector3 MinTsdFit();
        void TToFSDFit();
        void TToFT0Fit();
        void LikeliFit();
        void SelectHits(); // Not implemented yet
        void ComputeTSD(const TVector3&);
        void SetParLimits(const double*, const double dv=100., const double dt=25.);
        void Clear();
        void PrintVertices();
        void PrintVertex(const std::string&, const int &);

        void InitPars();
        double myAsymmGaus(const int&, const double&);
        double myPolN(const int&, const double&);
        double myExpo(const int&, const double&);

        TGraph *fLikeli;
        bool fUseGraph;
        int fNPMTs;
        std::vector<TVector3> fPos;
        std::vector<TVector3> fDir;
        HitCluster* fCl;
        HitCluster* fClSel;

        float fTSD;
        float fAveTTOF;

        double fTOFFSET;
        vector<vector<double>> fLim;
        vector<vector<double>> fVertex;

        double fLim_x[2];
        double fLim_y[2];
        double fLim_z[2];
        double fLim_t[2];
        

        std::vector<std::vector<double>> fPar;

        int    fMinPoint[3];
};
