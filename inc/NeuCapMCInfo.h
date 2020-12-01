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
    int grandID         =-1;
    int grandPDG        =-1;
    int ancestorPDG     =-1;
    bool isPriNeu       =false;
    float captureTime   =-99999.;
    float capturePos[3] ={-99999.};
    int nGammas         =0;
    float eGamma[20]    ={0.};
    float eGammaTotal   =0.;
    bool    isFound     =false;
};
