#include "VertexFit.h"

const float CWATER=21.5833;

// Long tank geometry
//const float TANK_R=371.0; 
//const float TANK_Y=1042.0/2.;

// Shor tank geometry
const float TANK_R=400.;
const float TANK_Y=300.;

const float EPSILON=1.e-5; // Just avoid precision issue

VertexFit::VertexFit() :
fTOFFSET(0.),
fUseGraph(false)
{
    fCl=NULL;
    fClSel=NULL;
    fPos.clear();
    fDir.clear();
    fVertex.clear();
    fLim.clear();
    for(int i=0; i<kNVtx; i++)
    {
        vector<double> v(kNAxis, 0.);
        fVertex.push_back( v );
    }

    for(int i=0; i<kNAxis; i++)
    {
        vector<double> v(kNEdge, 0.);
        fLim.push_back( v );
    }

    InitPars();
}

VertexFit::~VertexFit()
{
    if( !fCl ){ delete fCl; fCl=NULL; }
    if( !fClSel ){ delete fClSel; fClSel=NULL; }
    if( !fLikeli ){ delete fLikeli; fLikeli=NULL; }

    vector<vector<double>>().swap(fLim);
    vector<vector<double>>().swap( fVertex );
    vector<TVector3>().swap( fPos );
    vector<TVector3>().swap( fDir );
}


void VertexFit::LoadTemplate(std::string likelifile)
{
    TFile *f=TFile::Open(likelifile.c_str());  
    fLikeli=(TGraph*)f->Get("glikeli")->Clone();
    f->Close();
    fUseGraph=true;
}

void VertexFit::FitVertex(HitCluster *aCl)
{
    fCl=aCl; 
    this->Clear();

    //// Finding seed by fixed grid seach (very coase)
    TVector3 pre_pos=this->MinTsdFit();
    for(int i=eX; i<=eZ; i++){ fVertex[eSeed][i]=pre_pos(i); }
    fVertex[eSeed][eT]=this->GetAveTTOF();

    //// Using the result of the grid search, perform finner fit
    this->TToFSDFit();
    this->TToFT0Fit();

    // Selecting hits used for the following precise fit
    // (NOT IMPLEMENTED YET)
    ////this->SelectHits();

    // Apply precise fit
    this->LikeliFit();
    //this->PrintVertices();
}

void VertexFit::TToFSDFit()
{
    double pre[4]={0.}; 
    for(int i=eX; i<=eT; i++){ pre[i]=fVertex[eSeed][i]; }

    ResTFcn Fcn(this);
    Fcn.SetFlagTToFSDFit( true );
    MnUserParameters upar;
    for(int i=eX; i<=eZ; i++){ upar.Add(kSAxis[i], pre[i], 0.1); }
    this->SetParLimits(pre);
    for(int i=eX; i<=eZ; i++){ upar.SetLimits(kSAxis[i], fLim[i][eLow], fLim[i][eUp]); }

    extern int gErrorIgnoreLevel; gErrorIgnoreLevel = 1001;
    MnMigrad migrad(Fcn, upar);
    FunctionMinimum min = migrad();
    MnUserParameters res=min.UserParameters();
    for(int i=eX; i<=eZ; i++){ fVertex[eTToFSD][i]=res.Value(i); }

    int nHits=fCl->N();
    std::vector<float> tres;
    tres.reserve( nHits );
    int id=-9999999;
    double ave_ttof=0.;
    double dist=0.;
    for(int i=0; i<nHits; i++)
    {
        id=fCl->C(i)-1;
        for(int j=eX; j<=eZ; j++)
        {
            dist+=(fVertex[eTToFSD][j] - fPos[id](j))*(fVertex[eTToFSD][j] - fPos[id](j));
        }
        dist=TMath::Sqrt( dist );
        tres[i]=fCl->T(i) - dist/CWATER;
        ave_ttof+=tres[i];
    }
    ave_ttof=(double)ave_ttof/nHits;
    fVertex[eTToFSD][eT]=ave_ttof ;
}


void VertexFit::TToFT0Fit()
{
    double pre[4]={0.}; 
    for(int i=eX; i<=eT; i++){ pre[i]=fVertex[eTToFSD][i]; }

    ResTFcn Fcn(this);
    Fcn.SetFlagTToFT0Fit( true );
    MnUserParameters upar;
    for(int i=eX; i<=eT; i++){ upar.Add(kSAxis[i], pre[i], 0.1); }
    this->SetParLimits(pre);
    for(int i=eX; i<=eT; i++){ upar.SetLimits(kSAxis[i], fLim[i][eLow], fLim[i][eUp]); }

    extern int gErrorIgnoreLevel; gErrorIgnoreLevel = 1001;
    MnMigrad migrad(Fcn, upar);
    FunctionMinimum min = migrad();
    MnUserParameters res=min.UserParameters();
    for(int i=eX; i<=eT; i++){ fVertex[eTToFT0][i]=res.Value(i); }
}

void VertexFit::LikeliFit()
{
    double pre[4]={0.}; 
    for(int i=eX; i<=eT; i++){ pre[i]=fVertex[eTToFSD][i]; }

    ResTFcn Fcn(this);
    Fcn.SetFlagLikeliFit( true );
    MnUserParameters res;
    MnUserParameters upar;

    for(int i=eX; i<=eT; i++){ upar.Add(kSAxis[i], pre[i], 0.1); }
    this->SetParLimits(pre);
    for(int i=eX; i<=eZ; i++){ upar.SetLimits(kSAxis[i], fLim[i][eLow], fLim[i][eUp]); }

    // Suppressing print level of MINUIT2
    extern int gErrorIgnoreLevel; gErrorIgnoreLevel = 1001;

    // Should make them overwritable later
    const int nScan=3;
    const double tStep[nScan]      ={1., 0.2, 0.025};   // step size (ns)
    const double tEdgeL[nScan]     ={-10., -1.5, -0.4};    // scan range (ns)
    const double tEdgeU[nScan]     ={50., 1.5, 0.4};    // scan range (ns)
    const int    thrNpgrad[nScan]  ={4, 5, 4};
    int    nSteps[nScan]           ={0};
    for(int i=0; i<nScan; i++){ nSteps[i]=int( (tEdgeU[i]-tEdgeL[i])/tStep[i] ); }

    double minVtx[nScan][4]={0.}; 
    double mint0[nScan]={0.};
    double minFval[nScan]={0.};
    for(int iScan=0; iScan<nScan; iScan++){ minFval[iScan]=999999999.; }

    double t0=0.;
    double curFval=0.;
    double preFval=0.;
    int Npgrad=0;

    for(int iScan=0; iScan<nScan; iScan++)
    {
        if( iScan==0 ){ t0=pre[3]; }
        else          { t0=mint0[iScan-1]; }
        t0+=tEdgeL[iScan];
        Npgrad=0;

        for(int jp=0; jp<nSteps[iScan]; jp++)
        {
            if( jp!=0 ){ preFval=curFval; }

            upar.SetValue("t", t0);
            upar.Fix("t");

            MnMigrad migrad(Fcn, upar);
            FunctionMinimum min = migrad();
            curFval=min.Fval();
            if( curFval<minFval[iScan] )
            {
                minFval[iScan]  =curFval;
                mint0[iScan]    =t0;
                Npgrad          =0;
                res             =min.UserParameters();
                for(int k=0; k<4; k++){ minVtx[iScan][k]=res.Value(k); }
                fMinPoint[iScan]=jp;
            }

            if( jp!=0 && curFval-preFval>0. ){ Npgrad+=1; }
            if( Npgrad>=thrNpgrad[iScan] )
            { 
                //std::cout<<" [INFO] VertexFit::LikeliFit - Breaking" << iScan+1 
                //         <<" th loop at :" << jp <<std::endl;
                break;
            }
            t0+=tStep[iScan];
        }
    }
    for(int i=eX; i<=eT; i++)
    {
        fVertex[eNLLike][i]=minVtx[nScan-1][i];
    }
}


TVector3 VertexFit::MinTsdFit()
{
    bool  do_shrink=false;
    float Inc=1.; // Get smaller if do_shrink is true

    float max_scan_x=0.;
    float max_scan_y=0.;
    float max_scan_z=0.;

    TVector3 vtxTry(0., 0., 0.); // Current trial vertex
    TVector3 vtxMin(0., 0., 0.); // Vertex having minimum tsd
    TVector3 vtxPre(0., 0., 0.); // !!!! initial value is important
    
    float tryR2 =0.;
    float tsd   =0.;
    float tsdMin=9999999999.; // !!!! Initial value is important

    const float DWALL=25.;
    const float MAX_SCAN_XZ=2.*(TANK_R - DWALL);
    const float MAX_SCAN_Y =2.*(TANK_Y - DWALL);
    const int NDIVISIONS=10;
    const float STEP_Y =(float)MAX_SCAN_Y/NDIVISIONS;
    const float STEP_XZ=(float)MAX_SCAN_XZ/NDIVISIONS;

    while( Inc>0.01 )
    {
        max_scan_y=Inc*TANK_Y + vtxPre(1);
        max_scan_x=Inc*TANK_R + vtxPre(0);
        max_scan_z=Inc*TANK_R + vtxPre(2);

        if( max_scan_y>TANK_Y ){ max_scan_y=TANK_Y; }
        if( max_scan_x>TANK_R ){ max_scan_x=TANK_R; }
        if( max_scan_z>TANK_R ){ max_scan_z=TANK_R; }

        vtxTry(1)=-max_scan_y + vtxPre(1);
        while( vtxTry(1)<=max_scan_y )
        {
            vtxTry(0)=-max_scan_x + vtxPre(0);
            while( vtxTry(0)<=max_scan_x )
            {
                vtxTry(2)=-max_scan_z + vtxPre(2);
                while( vtxTry(2)<=max_scan_z )
                {
                    tryR2=vtxTry(0)*vtxTry(0) + vtxTry(2)*vtxTry(2);
                    if( TMath::Sqrt(tryR2)<TANK_R )
                    {
                        // Compute standard deviation of t-tof using
                        // all the hits of fCl (i.e. current hit cluster)
                        this->ComputeTSD(vtxTry);
                        tsd=this->GetTSD();

                        if( tsd<tsdMin )
                        {
                            tsdMin=tsd;
                            vtxMin=vtxTry;
                        }
                    }
                    vtxTry(2)+=Inc*STEP_XZ;
                }
                vtxTry(0)+=Inc*STEP_XZ;
            }
            vtxTry(1)+=Inc*STEP_Y;
        }

        if( do_shrink )
        {
            vtxPre=vtxTry;
            Inc/=2.;
        }
        else
        {
            break;
        }
    }
    return vtxMin;
}

void VertexFit::ComputeTSD(const TVector3 &v)
{
    float ave_ttof=0.;
    float tsd=0.;
    float dist=0.;
    int id=-9999;
    int nHits=fCl->N();
    std::vector<float> tres;
    tres.reserve( nHits );

    for(int i=0; i<nHits; i++)
    {
        id=fCl->C(i)-1;
        for(int j=0; j<3; j++)
        {
            dist+=(v(j) - fPos[id](j))*(v(j) - fPos[id](j));
        }
        dist=TMath::Sqrt( dist );
        tres[i]=fCl->T(i) - dist/CWATER;
        ave_ttof+=tres[i];
    }
    ave_ttof=(float)ave_ttof/nHits;
    
    for(int i=0; i<nHits; i++)
    {
        tsd+=(tres[i]-ave_ttof)*(tres[i]-ave_ttof);
    }
    tsd=(float)tsd/nHits;
    tsd=TMath::Sqrt( tsd );
   
    // Update private members 
    fTSD        =tsd;
    fAveTTOF    =ave_ttof;
}

double VertexFit::GetNLLike(const double *vtx)
{
    int nHits=fCl->N();
    double dist=0.;
    int id=-9999;
    double t_tof_t0=0.;
    double likelihood=0.;
    double nll=0.;

    for(int i=0; i<nHits; i++)
    {
        id=fCl->C(i)-1;
        for(int j=0; j<3; j++)
        {
            dist+=(vtx[j] - fPos[id](j))*(vtx[j] - fPos[id](j));
        }
        dist=TMath::Sqrt( dist );

        t_tof_t0=fCl->T(i) - dist/CWATER - vtx[3];
        if( t_tof_t0<-50. || t_tof_t0>250. )
        {
            likelihood=1.e-10;
        }
        else
        {
            if( fUseGraph ){ likelihood=fLikeli->Eval( t_tof_t0+10. ); }
            else           { likelihood=this->GetLogLikeli(t_tof_t0+10.); }

            if( likelihood<1.e-10 )
            {
                likelihood=1.e-10;
            }
        }
        nll-=TMath::Log( likelihood );
    }
    return nll;
}

double VertexFit::GetTToFT0(const double *vtx)
{
    int nHits=fCl->N();
    double dist=0.;
    int id=-9999;
    double t_tof_t0=0.;
    double gdn=0.;

    for(int i=0; i<nHits; i++)
    {
        id=fCl->C(i)-1;
        for(int j=0; j<3; j++)
        {
            dist+=(vtx[j] - fPos[id](j))*(vtx[j] - fPos[id](j));
        }
        dist=TMath::Sqrt( dist );

        t_tof_t0=fCl->T(i) - dist/CWATER - vtx[3];
        t_tof_t0/=0.9; 
        gdn+=TMath::Exp( -0.5*t_tof_t0*t_tof_t0 );
    }
    return gdn;
}

double VertexFit::GetTToFSD(const double *vtx)
{
    this->ComputeTSD(TVector3(vtx[0], vtx[1], vtx[2]));
    return (double)this->GetTSD();
}


void VertexFit::SetTrueVertex(const double *v)
{
    for(int i=eX; i<=eT; i++){ fVertex[eTrue][i]=v[i]; }
}

void VertexFit::SetParLimits(const double *v, const double dv, const double dt)
{
    const double fac=TANK_Y/TANK_R;
    fLim[eX][eLow]=v[eX]-dv;         fLim[eX][eUp]=v[eX]+dv;
    fLim[eY][eLow]=v[eY]-(dv*fac);   fLim[eY][eUp]=v[eY]+(dv*fac);
    fLim[eZ][eLow]=v[eZ]-dv;         fLim[eZ][eUp]=v[eZ]+dv;
    fLim[eT][eLow]=v[eT]-dt;         fLim[eT][eUp]=v[eT]+2.*dt;

    if( fLim[eX][eLow]<-TANK_R ){ fLim[eX][eLow]=-TANK_R; } 
    if( fLim[eY][eLow]<-TANK_Y ){ fLim[eY][eLow]=-TANK_Y; } 
    if( fLim[eZ][eLow]<-TANK_R ){ fLim[eZ][eLow]=-TANK_R; } 
    if( fLim[eX][eUp]>TANK_R ){ fLim[eX][eUp]=TANK_R; } 
    if( fLim[eY][eUp]>TANK_Y ){ fLim[eY][eUp]=TANK_Y; } 
    if( fLim[eZ][eUp]>TANK_R ){ fLim[eZ][eUp]=TANK_R; } 
}

void VertexFit::PrintVertex(const std::string &sName, const int &i)
{
    std::cout<<"  ----- " << sName 
             <<" : ("     << fVertex[i][eX]
             <<", "       << fVertex[i][eY]
             <<", "       << fVertex[i][eZ]
             <<", "       << fVertex[i][eT]
             <<")"        <<std::endl;

}
void VertexFit::PrintVertices()
{
    PrintVertex(std::string("True"), eTrue);
    PrintVertex(std::string("Seed"), eSeed);
    PrintVertex(std::string("TToFSD"), eTToFSD);
    PrintVertex(std::string("TToFT0"), eTToFT0);
    PrintVertex(std::string("NLLike"), eNLLike);
}

void VertexFit::Clear()
{
    for(unsigned int i=0; i<fVertex.size(); i++)
    {
        for(unsigned int j=0; j<fVertex[i].size(); j++)
        {
            fVertex[i][j]=0.;
        }
    }
}

void VertexFit::InitPars()
{
    fPar.clear();
    std::vector<double> Par1(4, 0.);
    Par1[0]=0.0260887;
    Par1[1]=-0.150878;
    Par1[2]=0.987932;
    Par1[3]=1.015;
    
    std::vector<double> Par2(8, 0.);
    Par2[0]=1.5;
    Par2[1]=0.00695033;
    Par2[2]=-0.00666912;
    Par2[3]=0.00322981;
    Par2[4]=-0.000850271;
    Par2[5]=0.000124309;
    Par2[6]=-9.48304e-06;
    Par2[7]=2.94241e-07;
    
    std::vector<double> Par3(7, 0.);
    Par3[0]=9;
    Par3[1]=0.000554756;
    Par3[2]=-4.31405e-05;
    Par3[3]=-2.25578e-07;
    Par3[4]=1.1608e-06;
    Par3[5]=-1.99993e-07;
    Par3[6]=1.13596e-08;
    
    std::vector<double> Par4(3, 0.);
    Par4[0]=16;
    Par4[1]=-7.95583;
    Par4[2]=-0.0471305;
    
    std::vector<double> Par5(7, 0.);
    Par5[0]=40;
    Par5[1]=0.000113132;
    Par5[2]=-5.64682e-06;
    Par5[3]=9.04757e-08;
    Par5[4]=4.17006e-10;
    Par5[5]=-9.20284e-12;
    Par5[6]=-3.13147e-13;

    fPar.push_back( Par1 );
    fPar.push_back( Par2 );
    fPar.push_back( Par3 );
    fPar.push_back( Par4 );
    fPar.push_back( Par5 );
}

double VertexFit::GetLogLikeli(const double &x)
{
    double likeli=0.1e-4;;
    if( x>-4.1 && x<1.5 ) // Asymmetric gaussian
    {
        likeli=myAsymmGaus(0, x);
    }
    else if( x>=1.5 && x<9. ) // 6th of polynomial
    {
        likeli=myPolN(1, x);
    }
    else if( x>=9. && x<16. )
    {
        likeli=myPolN(2, x);
    }
    else if( x>=16. && x<=40. )
    {
        likeli=myExpo(3, x);
    }
    else if( x>=40 && x<70. )
    {
    }
   return TMath::Log( likeli ); 
   //return likeli; 
}

double VertexFit::myAsymmGaus(const int &idx, const double &x)
{
	Double_t X=x-fPar[idx][1];
	Double_t sigma=X<0 ? fPar[idx][2] : fPar[idx][3];
    X=X/sigma;
	return fPar[idx][0]*TMath::Exp( -(X*X)/2.0 );
}

double VertexFit::myPolN(const int &idx, const double &x)
{
    double X=x-fPar[idx][0]; 
    double f=fPar[idx][1];
    for(unsigned int i=2; i<fPar[idx].size(); i++)
    {
        double xEn=1.;
        for(unsigned int j=1; j<=i-1; j++){ xEn*=X; }
        f+=fPar[idx][i]*xEn;
    }
    return f;
}

double VertexFit::myExpo(const int &idx, const double &x)
{
    double X=x - fPar[idx][0];
    return TMath::Exp(fPar[idx][1] + fPar[idx][2]*X);
}
