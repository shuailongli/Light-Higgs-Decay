#ifndef Model_H
#define Model_H

#include<iostream>
#include <cmath>
#include <complex>

#include "Parameters.h"

using namespace std;

class BSModel
{
private:
    
    //model type:   0, setting all couplings manuallys;
    //              1, 2HDM Type-I;
    //              2, 2HDM Type-II;
    //              3, 2HDM Type-L;
    //              4, 2HDM Type-F;
    //              5, MSSM;
    //              6, NMSSM;
    int model=6;
    
    //all masses are in the units of GeV
    //2HDM parameters
    double mA=0.5;          //CP-odd Higgs mass
    double mH=0.5;          //CP-even Higgs mass
    double mHC=600.;        //charged Higgs mass
    double tanb=10;
    double cosba=0;
    double lambdav=0;       //sqrt(lambda v^2)
    
    
    //normalized A-fermion-fermion couplings
    double xiAuu=1/tanb;    //normalized Auu couplings
    double xiAcc=1/tanb;    //normalized Acc couplings
    double xiAtt=1/tanb;    //normalized Att couplings
    
    double xiAdd=tanb;      //normalized Add couplings
    double xiAss=tanb;      //normalized Ass couplings
    double xiAbb=tanb;      //normalized Abb couplings
    
    double xiAee=tanb;      //normazlied Aee couplings
    double xiAmumu=tanb;    //normazlied Amumu couplings
    double xiAtata=tanb;    //normazlied Atatas couplings
    
    
    //normalized H-fermion-fermion couplings
    double xiHuu=cosba-sqrt(1-cosba*cosba)/tanb;    //normalized Huu couplings
    double xiHcc=cosba-sqrt(1-cosba*cosba)/tanb;    //normalized Hcc couplings
    double xiHtt=cosba-sqrt(1-cosba*cosba)/tanb;    //normalized Htt couplings
    
    double xiHdd=cosba+sqrt(1-cosba*cosba)*tanb;    //normalized Hdd couplings
    double xiHss=cosba+sqrt(1-cosba*cosba)*tanb;    //normalized Hss couplings
    double xiHbb=cosba+sqrt(1-cosba*cosba)*tanb;    //normalized Hbb couplings
    
    double xiHee=cosba+sqrt(1-cosba*cosba)*tanb;    //normazlied Hee couplings
    double xiHmumu=cosba+sqrt(1-cosba*cosba)*tanb;  //normazlied Hmumu couplings
    double xiHtata=cosba+sqrt(1-cosba*cosba)*tanb;  //normazlied Htatas couplings
    
    //H-V-V couplings for V= Z, W
    double xiHZZ=cosba, xiHWW=cosba;
    
    //H-HC-HC couplings, in the unit of 1/GeV
    double gHCC=-1/vev*(-lambdav*lambdav*(1/tanb-tanb)*sqrt(1-cosba*cosba)+ (2*mHC*mHC-mH*mH+2*lambdav*lambdav)*cosba);
    
    //MSSM parameters, mass in GeV
    double WinoMass=500, mueff=500;
    int Qchargino=-1;   //chargino charge in MSSM/NMSSM
    
    double lambdaSHH=0.3;   //Singlet-Hu-Hd coupling, dimensionless
    double P11=0.03;        //cos theta_P, CP-odd singlet and Higgs mixing
    double P12=-sqrt(1-0.03*0.03);  //-sin theta_P, =-sqrt(1-P11*P11)
    
    
    //Chargino Mass, m_chi1, m_chi2. Mass in GeV
    double CharginoMass[2]={564.476, 440.643};
    //angles of U & V that diagonalize Chargino mass matrix by UMV^T
    double thetaU=0.734911;
    double thetaV=0.835885;
    
    //Higgs effective couplings
    complex <double> Cgluon;    //A-gluon-gluon effective couplings, see (3) in 1612.06538
    complex<double> Cgamma;     //A-gamma-gamma effective couplings, see (3) in 1612.06538
    
    
    //effective Yukawa couplings in spectator model
    double YuAu;
    double YuAd;
    double YuAs;
    
    
    //A-meson mixing:   A-A, A-pi0, A-eta, A-eat'
    double mixAA;
    double mixA3;
    double mixAeta;
    double mixAetap;
    
    
    //U(3) generators in meson basis: pi+, pi-, pi0, K+, K-, K0, BK0, eta, eta'
    double Meson[9][3][3]={
        {{0, 1, 0}, {0, 0, 0}, {0, 0, 0}},
        {{0, 0, 0}, {1, 0, 0}, {0, 0, 0}},
        {{1/sqrt(2.), 0, 0}, {0, -1/sqrt(2.), 0}, {0, 0, 0}},
        {{0, 0, 1}, {0, 0, 0}, {0, 0, 0}},
        {{0, 0, 0}, {0, 0, 0}, {1, 0, 0}},
        {{0, 0, 0}, {0, 0, 1}, {0, 0, 0}},
        {{0, 0, 0}, {0, 0, 0}, {0, 1, 0}},
        {{1/sqrt(6.)*coseta-1/sqrt(3.)*sineta, 0, 0},
            {0, 1/sqrt(6.)*coseta-1/sqrt(3.)*sineta, 0},
            {0, 0, -2/sqrt(6.)*coseta-1/sqrt(3.)*sineta}},
        {{1/sqrt(6.)*sineta+1/sqrt(3.)*coseta, 0, 0},
            {0, 1/sqrt(6.)*sineta+1/sqrt(3.)*coseta, 0},
            {0, 0, -2/sqrt(6.)*sineta+1/sqrt(3.)*coseta}}};
    
    
    
    //A-meson-meson-meson couplings in Chiral Perturbation Theory
    //pi1, pi2, pi3: IDs for 9 mesons
    //                  1:  pi+
    //                  2:  pi-
    //                  3:  pi0
    //                  4:  K+
    //                  5:  K-
    //                  6:  K0
    //                  7:  BK0
    //                  8:  eta
    //                  9:  eta'
    double cupAPiPiPi(int pi1, int pi2, int pi3);
    
    //Meson-meson-meson-meson couplings in Chiral Perturbation Theory
    //pi1, pi2, pi3, pi4: IDs for 9 mesons
    //                  1:  pi+
    //                  2:  pi-
    //                  3:  pi0
    //                  4:  K+
    //                  5:  K-
    //                  6:  K0
    //                  7:  BK0
    //                  8:  eta
    //                  9:  eta'
    double cupPiPiPiPi(int pi1, int pi2, int pi3, int pi4);
    
    
    
    //call this function to update parameters whenever a change to the input parameters is made
    void Initialize();
    
    
    
    
public:
    
    BSModel(){
        Initialize();
    }
    
    virtual ~BSModel(){}
    
    
    //********************************************************************************//
    //********************************* Set up model *********************************//
    //********************************************************************************//
    
    void set_model(int model_in){
        if(model_in<0 || model_in>6){
            cout<<"*Warning: Model Type must be an integer between 0-6!"<<endl;
            exit(0);
        }
        else{
        model=model_in;
        Initialize();
        }
    }
    
    //input mA in GeV
    void set_mA(double mA_in){
        mA=mA_in;
        Initialize();
    }
    
    //input mH in GeV
    void set_mH(double mH_in){
        mH=mH_in;
        Initialize();
    }
    
    //input mHC in GeV
    void set_mHC(double mH_in){
        mHC=mH_in;
        Initialize();
    }
    
    //input lambdav in GeV
    void set_lambdav(double mH_in){
        lambdav=mH_in;
        Initialize();
    }
    
    void set_tanb(double tanb_in){
        tanb=tanb_in;
        Initialize();
    }
    
    void set_cosba(double cosba_in){
        cosba=cosba_in;
        Initialize();
    }
    
    void set_P11(double P11_in){
        P11=P11_in;
        Initialize();
    }
    
    //input WinoMass in GeV
    void set_WinoMass(double WinoMass_in){
        WinoMass=WinoMass_in;
        Initialize();
    }
    
    //input mueff in GeV
    void set_mueff(double mueff_in){
        mueff=mueff_in;
        Initialize();
    }
    
    void set_lambdaSHH(double lambdaSHH_in){
        lambdaSHH=lambdaSHH_in;
        Initialize();
    }
    
    //initialize xiAtt, xiAbb, xiAcc, xiAss, xiAuu, xiAdd, xiAtata, xiAmumu, xiAee
    //can only be invoked for model==0
    void set_AffCoupling(double,double,double,double,double,double,double,double,double);
    
    //initialize xiHtt, xiHbb, xiHcc, xiHss, xiHuu, xiHdd, xiHtata, xiHmumu, xiHee
    //can only be invoked for model==0
    void set_HffCoupling(double,double,double,double,double,double,double,double,double);
    
    //initialize xiHZZ, xiHWW
    //can only be invoked for model==0
    void set_HVVCoupling(double,double);
    
    //initialize gHCC
    //can only be invoked for model==0
    void set_gHCC(double);
    
    
    
    //********************************************************************************//
    //*************************** functions for CP-odd Higgs *************************//
    //********************************************************************************//
    
    //diphoton decay, in GeV
    double DecayWidth_Agaga();          // A-> gamma gamma decay width
    
    //leptonic decay, in GeV
    double DecayWidth_Aee();            // A-> e+ e- decay width
    double DecayWidth_Amumu();          // A-> mu+ mu- decay width
    double DecayWidth_Atautau();        // A-> tau+ tau- decay width
    
    //hadronic decay in Chiral limit, mA<1 GeV. in unit of GeV
    double DecayWidth_APiPiPi();        // A-> pi pi pi
    double DecayWidth_AEtaPiPi();       // A-> eta pi pi
    double DecayWidth_AEtapPiPi();      // A-> eta' pi pi
    double DecayWidth_AEtaEtaPi();      // A-> eta eta pi
    double DecayWidth_AKKPi();          // A-> K K pi
    double DecayWidth_AGammaPiPi();     // A-> gamma pi pi
    
    //hadronic decay in spectator model, 1 GeV<mA<3 GeV, in GeV
    double DecayWidth_AEtaEtapPi();     // A-> eta eta' pi
    double DecayWidth_AEtapEtapPi();    // A-> eta' eta' pi
    double DecayWidth_AEtaEtaEta();     // A-> eta eta eta
    double DecayWidth_AEtaEtaEtap();    // A-> eta eta eta'
    double DecayWidth_AEtaEtapEtap();   // A-> eta eta' eta'
    double DecayWidth_AEtapEtapEtap();  // A-> eta' eta' eta'
    double DecayWidth_AEtaKK();         // A-> eta K K
    double DecayWidth_AEtapKK();        // A-> eta' K K
    
    //partonic hardon decay, mA>3 GeV, in GeV
    double DecayWidth_AQuark();         // A-> quark quark
    double DecayWidth_AGluon();         // A-> gluon gluon
    
    
    //Total hadronic decay width, in GeV
    double DecayWidth_ATotalHadron();
    
    //Total leptonic decay width, in GeV
    double DecayWidth_ATotalLepton();
    
    //Total decay width, in GeV
    double DecayWidth_ATotal();
    
    
    //Decay Branching Fractions
    double BranchRatio_Agaga();         // A-> gamma gamma branching fraction
    
    //leptonic decay Branching Fractions
    double BranchRatio_Aee();           // A-> e+ e- decay branching fraction
    double BranchRatio_Amumu();         // A-> mu+ mu- decay branching fraction
    double BranchRatio_Atautau();       // A-> tau+ tau- decay branching fraction
    
    //Total hadronic decay Branching Fractions
    double BranchRatio_ATotalHadron();
    
    //partonic hardon decay, mA>3 GeV
    double BranchRatio_AQuark();        // A-> quark quark branching fraction
    double BranchRatio_AGluon();        // A-> gluon gluon branching fraction
    
    
    //********************************************************************************//
    //*************************** functions for CP-even Higgs ************************//
    //********************************************************************************//
    
    //diphoton decay, in GeV
    double DecayWidth_Hgaga();          // H-> gamma gamma decay width

    //leptonic decay, in GeV
    double DecayWidth_Hee();            // H-> e+ e- decay width
    double DecayWidth_Hmumu();          // H-> mu+ mu- decay width
    double DecayWidth_Htautau();        // H-> tau+ tau- decay width
    
    //meson decay, in GeV
    double DecayWidth_HPiPi();          // H-> pi pi
    double DecayWidth_HKK();            // H-> Kaon Kaon
    double DecayWidth_HPiPiPiPi();      // H-> pi pi pi pi
    
    //partonic hardon decay, in GeV
    double DecayWidth_Hcc();            // H-> charm charm
    double DecayWidth_Hss();            // H-> strange strange
    double DecayWidth_HGluon();         // H-> gluon gluon
    
    //Total decay width, in GeV
    double DecayWidth_HTotal();
    
    
    //Decay Branching Fractions
    double BranchRatio_Hgaga();         // H-> gamma gamma branching fraction
    
    //leptonic decay Branching Fractions
    double BranchRatio_Hee();           // H-> e+ e- decay branching fraction
    double BranchRatio_Hmumu();         // H-> mu+ mu- decay branching fraction
    double BranchRatio_Htautau();       // H-> tau+ tau- decay branching fraction
    
    //meson decay branching fractions
    double BranchRatio_HPiPi();         // H-> pi pi
    double BranchRatio_HKK();           // H-> Kaon Kaon
    double BranchRatio_HPiPiPiPi();     // H-> pi pi pi pi
    
    //partonic hardon decay, mA>3 GeV
    double BranchRatio_Hcc();           // H-> charm charm
    double BranchRatio_Hss();           // H-> strange strange
    double BranchRatio_HGluon();        // H-> gluon gluon branching fraction

};


#endif
