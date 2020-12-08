#pragma once

struct Track_t
{
    int thisID      =-1;
    int parentID    =-1;
    int thisPDG     =-1;
    float time      =-99999.;
    float posSrt[3] ={-999999.};
    float posEnd[3] ={-999999.};
    float energy    =-1.;
    float dirSrt[3] ={-99999.};
};

struct MCNeuCap_t 
{
    int nucleusPDG      =-1;
    int parentID        =-1;
//    int grandID         =-1;
//    int grandPDG        =-1;
//    int ancestorPDG     =-1;
//    bool isPriNeu       =false;
//
    float captureTime   =-99999.;
    float capturePos[3] ={-99999.};
    float captureDwall  =0.;
    int nGammas         =0;
    float eGamma[20]    ={0.};
    float eGammaTotal   =0.;

    int nAncestor       =0;
    float antrVtx_x[50] ={-9999.};
    float antrVtx_y[50] ={-9999.};
    float antrVtx_z[50] ={-9999.};
    float antrVtx_t[50] ={-9999.};
    float antrDir_x[50] ={-9999.};
    float antrDir_y[50] ={-9999.};
    float antrDir_z[50] ={-9999.};
    float antrEnergy[50]={-9999.};
    int antrPrntID[50]  ={-9999};
    int antrPDG[50]     ={-9999};
};

struct MichelE_t
{
    int grandID         =-1;
    int parentPDG       =-1;
    float Pos[3]        ={-9999.};
    float Time          =-9999.;
    float Dwall         =-9999.;
    float energy        =-9999.;
};
