#ifndef _DEFINITIONS_HH_
#define _DEFINITIONS_HH_

#include"libs.hh"
#include"R3BGladFieldMap.h"
#include"R3BTPropagator.h"
#include"R3BTrackingParticle.h"
#include <TMath.h>

//------------- Parameters of the training sample tracks -----------
const UInt_t    NEVENTS             = 6000; //Number of tracks for the training sample
const Double_t  AMU                 = 0.9314940038; //GeV 90 from AME 2016
const Double_t  BEAM_Z              = 50.0000000000; //primary Z - 124Sn
const Double_t  BEAM_A              = 124.000000000;//Primary A - 124Sn
const Double_t  BEAM_MASS_GEV_C2    = 123.905279 * AMU;// 124Sn mass in GeV/c2
const Double_t  BEAM_DTHETA_DEG     = 0.03*180./3.14;//max spread of theta angle in deg (+-)
const Double_t  BEAM_RADIUS_CM      = 3;//radius of the beamspot in cm 
const Double_t  BEAM_EKIN_GEV_U     = 0.900;//per unit in GeV - 120Sn
const Double_t  SCALE_MIN           = 0.15;// scaling for subtraction: p_min = Ptot-Ptot*scale_min
const Double_t  SCALE_MAX           = 0.15; // scaling for aadition:   p_max = Ptot+Ptot*scale_max
const Double_t  GLAD_CURRENT        = -2983.0; //run 490

//Nominal detector positions 06/11/2022 based on the data from Ivana and Eleonora
//Positions in [cm] w.r.t. Turning Point (TP) of GLAD along 18 deg line
const Double_t Angle        = 18.0 * TMath::Pi() / 180.;//Central turning angle of GLAD
const Double_t Angle_det    = TMath::Pi()/2. - Angle;//angle of the detector
const Double_t Target_to_GLAD_flange = 176.65;
const Double_t Target_to_MWPC = 164.85;
const Double_t Target_to_TP = Target_to_GLAD_flange + 164.;
const Double_t TP_to_F10 = 381.;
const Double_t TP_to_F11 = 1347.2;
const Double_t TP_to_F12 = TP_to_F11 + 20;
const Double_t TP_to_TOFD = TP_to_F12 + 14;

const UInt_t  FIT_VALUE = 0; //0=PoQ, 1=TX0, 2=TY0, 3=FlightPath, 4=TX1 after GLAD, 5=TY1 after glad, 6=Y after glad

//------------- Set Multi Dimensional Fit Class for Tracking -----------

const UInt_t   nVars   	= 8;
const UInt_t   nReduct 	= 0;//numbers of vars to ignore in PCA 
const UInt_t   Polynom_Type     = 0; //0=kChebyshev,1 =kLegendre(fails for MDFWrapper, don't use!),2-kMonomials 
const UInt_t   Power   	        = 8;
const UInt_t   MDFMaxFunctions  = 20000;
const UInt_t   MDFMaxStudy      = 400000;
const UInt_t   MDFMaxTerms      = 80; //PoQ
const UInt_t   MDFPowerLimit    = 1;
const Double_t MDFMinAngle      = 0.1;
const Double_t MDFMaxAngle      = 0.1;
const Double_t MDFMinError      = 1e-19;

//------------- Plotting ranges for the trajectories in [cm] -------------

const  Double_t Ymin = -100;
const  Double_t Ymax =  100;
const  Double_t Xmin = -1050;
const  Double_t Xmax =  1050;
const  Double_t Zmin = -300;
const  Double_t Zmax =  1800;
const  UInt_t NbinsX = 1000;
const  UInt_t NbinsY = 1000;
const  UInt_t NbinsZ = 1000;

//-------------- Data container for a given track --------------

struct Data
{
    Double_t edata[nVars];//data before PCA
    Double_t pdata[nVars - nReduct];//output data from PCA
    Double_t value;//main fit value
    R3BTrackingParticle primary;//original primary particle
};

//------------ Declaration of all fucntions ------------ 

R3BTrackingParticle Generate_Primary();
bool ExtractData(Data &d, R3BTPropagator* fPropagator);
void Print_MDF_params(TMultiDimFit* mdf);
void Save_MDF_params(TMultiDimFit* mdf, const char * filename);
void Save_PCA_params(TPrincipal* pca, const char * filename);
void run_MDF(Bool_t use_PCA, R3BTPropagator * fPropagator);

#endif




