#pragma once

namespace GL
{
    const int kPDG_ELECTRON         =11;
    const int kPDG_MuMinus          =13;
    const int kPDG_GAMMA            =22;
    const int kPDG_NEUTRON          =2112;
    const int kPDG_NUCLEUS_OFFSET   =(int)1.e9;
    const int kPDG_PROTON           =2212;
    const int kPDG_PiPlus           =211;
    const int kPDG_PiMinus          =-211;

    const int kMaxTrueNCap=500;
    const int kMaxRecoNCap=500;
    
    const int kMaxTrueMiE=50;
    const int kMaxTrueAntr=10;

    enum ENCapType{ eDirectPriNeu=0, 
                    eNeutronSecNeu, 
                    eProtonSecNeu,
                    ePiPlusSecNeu,
                    ePiMinusSecNeu,
                    eMuMinusCapNeu,
                    eOtherSecNeu
                  }; 

    const int kNoiseDN      =1;
    const int kNoiseAP      =2;
    const float TWID        =70.;
    const int   N10MAX      =70;
    const int   kNhitsMAX    =90;
    const int   N10THR      =6;
};
