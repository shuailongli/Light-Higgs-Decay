#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <math.h> 
#include <complex.h>
#include <cmath>

#include "Model.h"
#include "ToolFunctions.h"

/*
#include "./data_base/abs_Delta_K.h"
#include "./data_base/abs_Delta_Pi.h"
#include "./data_base/abs_Gamma_K.h"
#include "./data_base/abs_Gamma_Pi.h"
#include "./data_base/abs_Theta_K.h"
#include "./data_base/abs_Theta_Pi.h"

#include "./data_base/arg_Delta_K.h"
#include "./data_base/arg_Delta_Pi.h"
#include "./data_base/arg_Gamma_K.h"
#include "./data_base/arg_Gamma_Pi.h"
#include "./data_base/arg_Theta_K.h"
#include "./data_base/arg_Theta_Pi.h"

#include "./data_base/alphasQ.h"
*/
using namespace std;

complex <double> i(0, 1), comp;

complex<double> AQ12H(double mphi, double mi){
        double xi=mphi*mphi/(4*mi*mi);
        complex<double> fxi=0;
        if (xi>1)
        {
            fxi= -0.25*pow(log( (1+sqrt(1-1/xi)) / (1-sqrt(1-1/xi)) )-i*M_PI,2);
        }
        if (xi<=1)
        {
            fxi=pow(asin(sqrt(xi)),2);
        }
        complex<double> aqr=(xi+(xi-1)*fxi)/xi/xi;// there is no factor 2 here!!!
        return aqr;
    }

complex<double> AQ1H(double mphi, double mi){
        double xi=mphi*mphi/(4*mi*mi);
        complex<double> fxi=0;
        if (xi>1)
        {
            fxi= -0.25*pow(log( (1+sqrt(1-1/xi)) / (1-sqrt(1-1/xi)) )-i*M_PI,2);
        }
        if (xi<=1)
        {
            fxi=pow(asin(sqrt(xi)),2);
        }
        complex<double> aqr=-0.5*(2*xi*xi+3*xi+3*(2*xi-1)*fxi)/xi/xi;// there is no factor 2 here!!!
        return aqr;
    }

 
complex<double> AQ0H(double mphi, double mi){
        double xi=mphi*mphi/(4*mi*mi);
        complex<double> fxi=0;
        if (xi>1)
        {
            fxi= -0.25*pow(log( (1+sqrt(1-1/xi)) / (1-sqrt(1-1/xi)) )-i*M_PI,2);
        }
        if (xi<=1)
        {
            fxi=pow(asin(sqrt(xi)),2);
        }
        complex<double> aqr=-0.5*(xi-fxi)/xi/xi;// there is no factor 2 here!!!
        return aqr;
    }   
double Gamma_pipi(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie){
    complex <double> i(0, 1), comp;
    double pi=3.1415926535;
    double mpion=0.13497;//GeV
    //double vev=246.2;//GeV
    double fcon=3/vev/vev/32.0/pi;

    //double mu=0.0022;
    //double md=0.0047;
    //double mc=1.28;
    //double ms=0.096;
    //double mb=4.18;
    //double mt=173.1;

    complex<double> xig=3./2.*(xic*AQ12H(coe,mc)+xib*AQ12H(coe,mb)+xit*AQ12H(coe,mt));

    // cout<<xig<<" ";


    
    double beta_pi=sqrt(1-4*mpion*mpion/coe/coe);
    // Gamma_Pi=(abs_Gamma_Pi(coe),0i)*exp(i*arg_Gamma_Pi(coe));
    // Delta_Pi=(abs_Delta_Pi(coe),0i)*exp(i*arg_Delta_Pi(coe));
    // Theta_Pi=(abs_Delta_Pi(coe),0i)*exp(i*arg_Theta_Pi(coe));
    complex<double> Gamma_Pi=abs_Gamma_Pi(coe)*exp(i*arg_Gamma_Pi(coe));
    complex<double> Delta_Pi=abs_Delta_Pi(coe)*exp(i*arg_Delta_Pi(coe));
    complex<double> Theta_Pi=abs_Theta_Pi(coe)*exp(i*arg_Theta_Pi(coe));

    double GPP=1e-40;
    
    GPP=fcon/coe*beta_pi*pow(abs(Gamma_Pi *(mu*xiu+md*xid)/(mu+md) + xis*Delta_Pi + xig*(Theta_Pi-Delta_Pi-Gamma_Pi)*2.0/27.0),2);
    if (coe<2*mpion){
        GPP=1e-40;
    }
    if (coe>2){
        GPP=1e-40;
    }  
    return GPP;
    
}

// KK
double Gamma_KK(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie){
    complex <double> i(0, 1), comp;
    double pi=3.1415926535;
    double mkaon=0.498;//GeV
    //double vev=246.2;//GeV
    double fcon=3/vev/vev/32.0/pi;
    //double mu=0.0022;
    //double md=0.0047;
    //double mc=1.28;
    //double ms=0.096;
    //double mb=4.18;
    //double mt=173.1;

    complex<double> xig=3.0/2.*(xic*AQ12H(coe,mc)+xib*AQ12H(coe,mb)+xit*AQ12H(coe,mt));


    complex<double>  Gamma_K,Delta_K,Theta_K;

    double beta_k=sqrt(1-4*mkaon*mkaon/coe/coe);
    // Gamma_Pi=(abs_Gamma_Pi(coe),0i)*exp(i*arg_Gamma_Pi(coe));
    // Delta_Pi=(abs_Delta_Pi(coe),0i)*exp(i*arg_Delta_Pi(coe));
    // Theta_Pi=(abs_Delta_Pi(coe),0i)*exp(i*arg_Theta_Pi(coe));
    Gamma_K=abs_Gamma_K(coe)*exp(i*arg_Gamma_K(coe));
    Delta_K=abs_Delta_K(coe)*exp(i*arg_Delta_K(coe));
    Theta_K=abs_Theta_K(coe)*exp(i*arg_Theta_K(coe));
    
    double GKK=1e-40;
    GKK=1/vev/vev/8.0/pi/coe*beta_k*pow(abs(Gamma_K * (mu*xiu+md*xid)/(mu+md)+ xis*Delta_K  + xig*(Theta_K-Gamma_K -Delta_K)*2.0/27.0),2);
    if (coe<2*mkaon)
    {
        GKK=1e-40;
    }
    if (coe>2)
    {
        GKK=1e-40;
    }
    return GKK;
}

// 4pi
double Gamma_4pi(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie,  double c4pi){
    double mpion=0.13497;//GeV
    double beta_4pi=sqrt(1-4*2*mpion*2*mpion/coe/coe);   
    double G4pi=0;
    //double mu=0.0022;
    //double md=0.0047;
    //double mc=1.28;
    //double ms=0.096;
    //double mb=4.18;
    //double mt=173.1;
    double xigg=abs(xic*AQ12H(coe,mc)+xib*AQ12H(coe,mb)+xit*AQ12H(coe,mt)+xis*AQ12H(coe,ms)+xiu*AQ12H(coe,mu)+xid*AQ12H(coe,md))/2.0;
    // c4pi=5.1E-9;
    // cout<<xigg<<endl;


    G4pi=c4pi*beta_4pi*pow(coe,3)*xigg*xigg;
    if (coe<2*2*mpion)
    {
        G4pi=1e-40;
    }
    if (coe>2)
    {
        G4pi=1e-40;
    }
    return G4pi;
}

// ee
double Gamma_ee(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie){
    //double vev=246.2;//GeV
    //double me=0.511/1000.0;//GeV
    double pi=3.1415926535;
    double Gamma_e=1e-40;

    Gamma_e=me*me/(8*pi*vev*vev)*coe*pow(1-4*me*me/coe/coe,1.5)*pow(xie,2);
    if (coe<2*me)
    {
        Gamma_e=1e-40;
    }
    return Gamma_e;
}

// muon
double Gamma_mumu(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie){
    //double vev=246.2;//GeV
    //double mmu=0.1057;//GeV
    double pi=3.1415926535;
    double Gamma_mu=1e-40;

    Gamma_mu=mmu*mmu/(8*pi*vev*vev)*coe*pow(1-4*mmu*mmu/coe/coe,1.5)*pow(ximu,2);
    if (coe<2*mmu)
    {
        Gamma_mu=1e-40;
    }
    // cout<<ximu<<endl;
    return Gamma_mu;
}
// 

// tautau
double Gamma_tautau(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie){
    //double vev=246.2;//GeV
    double pi=3.1415926535;
    //double mtau=1.7768;
    double Gamma_tau=1e-40;

    Gamma_tau=mtau*mtau/(8*pi*vev*vev)*coe*pow(1-4*mtau*mtau/coe/coe,1.5)*pow(xitau,2);
    if (coe<2*mtau)
    {
        Gamma_tau=1e-40;
    }
    // cout<<Gamma_tau<<" ";
    return Gamma_tau;
    
}
//
// cc
double Gamma_cc(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie){
    //double vev=246.2;//GeV
    double pi=3.1415926535;
    //double mc=1.2751;//GeV
    double mD = 1.864843;//GeV
    double Gamma_c=1e-40;

    Gamma_c=3*mc*mc/(8*pi*vev*vev)*coe*pow(1-4*mD*mD/coe/coe,1.5)*pow(xic,2);
    if (coe<2*mD)
    {
        Gamma_c=1e-40;
    }
    return Gamma_c;
}
// 
// ss
double Gamma_ss(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie){
    //double vev=246.2;//GeV
    double pi=3.1415926535;
    //double ms=0.096;
    double mkaon=0.498;//GeV
    double Gamma_s=1e-40;

    Gamma_s=3*ms*ms/(8*pi*vev*vev)*coe*pow(1-4*mkaon*mkaon/coe/coe,1.5)*pow(xis,2);
        
    if (coe<2)
    {
        Gamma_s=1e-40;
    }
    return Gamma_s;
}

// gg
double Gamma_gg(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie){

    //double mc=1.28;//GeV
    // double mc=1.869;//GeV
    //double ms=0.096;
    //double mtau=1.7768;
    //double mt=173.2;
    //double mb= 4.18;
    //double mu=0.0022;
    //double md=0.0047;
    //double vev=246.2;//GeV
    double pi=3.1415926535;
    double Gamma_g=1e-40;
    if (coe<2)
    {
        return 1e-40;
    }
    // from 1809.01876

    Gamma_g=pow(alphasQ(coe),2)*pow(coe,3)/(32*pow(pi,3)*vev*vev)*pow(abs(xic*AQ12H(coe,mc)+xib*AQ12H(coe,mb)+xit*AQ12H(coe,mt)+xis*AQ12H(coe,ms)+xiu*AQ12H(coe,mu)+xid*AQ12H(coe,md)),2);
    return Gamma_g;
}



// gammagamma
double Gamma_gaga(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie,double xiz, double xiw, double mch,double HCC){

    //double mc=1.28;//GeV
    // double mc=1.869;//GeV
    //double ms=0.096;
    //double mtau=1.7768;
    //double mt=173.2;
    //double mb= 4.18;
    //double mu=0.0022;
    //double md=0.0047;
    //double vev=246.2;//GeV
    double pi=3.1415926535;
    double mw=80.4;
    double mz=91.187;
    double Gamma_ga=1e-40;
    
    // if (coe<2)
    // {
    //     return 1e-40;
    // }
    //
    complex <double> xiHgaga_q=4.0/3.0*xic*AQ12H(coe,mc)+1.0/3.0*xib*AQ12H(coe,mb)+4.0/3.0*xit*AQ12H(coe,mt)+1.0/3.0*xis*AQ12H(coe,ms)+4.0/3.0*xiu*AQ12H(coe,mu)+1.0/3.0*xid*AQ12H(coe,md) ;
    complex <double>xiHgaga_w=xiw*AQ1H(coe,mw);
    complex <double> xiHgaga_C= -vev/2.0/mch/mch*HCC*AQ0H(coe,mch)  ;


    Gamma_ga=pow(1.0/137.0,2)*pow(coe,3)/(64*pow(pi,3)*vev*vev)*pow(abs(xiHgaga_C+xiHgaga_q+xiHgaga_w),2);
    return Gamma_ga;
}

double llp_bsm_src(double coe, double xit, double xib, double xic, double xis, double xiu, double xid, double xitau, double ximu, double xie,double xiz, double xiw, double mch,double HCC, \
      double & GammaT, double & BrPP, double & BrKK, double & Br4pi, double & Brmu, double & Bre,double & Brtau, double & Brc, double & Brs, double & Brgg, double & Brgaga){
    complex <double> i(0, 1), comp;

    double c4pi=5.1E-9;
    // double coe=2;
    double GPP=Gamma_pipi(coe,xit,xib,xic,xis,xiu,xid,xitau,ximu,xie);
    double GKK=Gamma_KK(coe,xit,xib,xic,xis,xiu,xid,xitau,ximu,xie);
    double G4pi=Gamma_4pi(coe,xit,xib,xic,xis,xiu,xid,xitau,ximu,xie,c4pi);
    double Gmu=Gamma_mumu(coe,xit,xib,xic,xis,xiu,xid,xitau,ximu,xie);
    double Ge=Gamma_ee(coe,xit,xib,xic,xis,xiu,xid,xitau,ximu,xie);
    double Gtau=Gamma_tautau(coe,xit,xib,xic,xis,xiu,xid,xitau,ximu,xie);
    double Gc=Gamma_cc(coe,xit,xib,xic,xis,xiu,xid,xitau,ximu,xie);
    double Gs=Gamma_ss(coe,xit,xib,xic,xis,xiu,xid,xitau,ximu,xie);
    double Ggg=Gamma_gg(coe,xit,xib,xic,xis,xiu,xid,xitau,ximu,xie);
    double Ggaga=Gamma_gaga(coe,xit,xib,xic,xis,xiu,xid,xitau,ximu,xie, xiz,xiw, mch,HCC);

    double length;
    GammaT=GPP+GKK+G4pi+Gmu+Ge+Gtau+Gc+Gs+Ggg+Ggaga;
    length=6.58*pow(10,-25)*2.998*pow(10,8)/GammaT;
    BrPP=GPP/GammaT;
    BrKK=GKK/GammaT;
    Br4pi=G4pi/GammaT;
    Brmu=Gmu/GammaT;
    Bre=Ge/GammaT;
    Brtau=Gtau/GammaT;
    Brc=Gc/GammaT;
    Brs=Gs/GammaT;
    Brgg=Ggg/GammaT;
    Brgaga=Ggaga/GammaT;
    // cout<<setw(5)<<coe<<" "<<setw(10)<<GPP;
    // cout<<" "<<setw(10)<<GKK;
    // cout<<" "<<setw(10)<<G4pi;
    // cout<<" "<<setw(10)<<Gmu;
    // cout<<" "<<setw(10)<<Gtau;
    // cout<<" "<<setw(10)<<Gc;
    // cout<<" "<<setw(10)<<Gs;
    // cout<<" "<<setw(10)<<Ggg;
    // cout<<" "<<setw(10)<<Ge;
    // cout<<" "<<setw(10)<<GammaT;
    // cout<<" "<<setw(10)<<length;
    // cout<<endl;
    return length;
}



// H-> gamma gamma decay width
double BSModel::DecayWidth_Hgaga(){
    return Gamma_gaga(mH,xiHtt,xiHbb,xiHcc,xiHss,xiHuu,xiHdd,xiHtata,xiHmumu,xiHee, xiHZZ,xiHWW, mHC,gHCC);
}

// H-> e+ e- decay width
double BSModel::DecayWidth_Hee(){
    return Gamma_ee(mH,xiHtt,xiHbb,xiHcc,xiHss,xiHuu,xiHdd,xiHtata,xiHmumu,xiHee);
}

// H-> mu+ mu- decay width
double BSModel::DecayWidth_Hmumu(){
    return Gamma_mumu(mH,xiHtt,xiHbb,xiHcc,xiHss,xiHuu,xiHdd,xiHtata,xiHmumu,xiHee);
}

// H-> tau+ tau- decay width
double BSModel::DecayWidth_Htautau(){
    return Gamma_tautau(mH,xiHtt,xiHbb,xiHcc,xiHss,xiHuu,xiHdd,xiHtata,xiHmumu,xiHee);
}

// H-> pi pi
double BSModel::DecayWidth_HPiPi(){
    return Gamma_pipi(mH,xiHtt,xiHbb,xiHcc,xiHss,xiHuu,xiHdd,xiHtata,xiHmumu,xiHee);
}

// H-> Kaon Kaon
double BSModel::DecayWidth_HKK(){
    return Gamma_KK(mH,xiHtt,xiHbb,xiHcc,xiHss,xiHuu,xiHdd,xiHtata,xiHmumu,xiHee);
}

// H-> pi pi pi pi
double BSModel::DecayWidth_HPiPiPiPi(){
    double c4pi=5.1E-9;
    return Gamma_4pi(mH,xiHtt,xiHbb,xiHcc,xiHss,xiHuu,xiHdd,xiHtata,xiHmumu,xiHee,c4pi);
}

// H-> charm charm
double BSModel::DecayWidth_Hcc(){
    return Gamma_cc(mH,xiHtt,xiHbb,xiHcc,xiHss,xiHuu,xiHdd,xiHtata,xiHmumu,xiHee);
}

// H-> strange strange
double BSModel::DecayWidth_Hss(){
    return Gamma_ss(mH,xiHtt,xiHbb,xiHcc,xiHss,xiHuu,xiHdd,xiHtata,xiHmumu,xiHee);
}

// H-> gluon gluon
double BSModel::DecayWidth_HGluon(){
    return Gamma_gg(mH,xiHtt,xiHbb,xiHcc,xiHss,xiHuu,xiHdd,xiHtata,xiHmumu,xiHee);
}

//Total decay width, in GeV
double BSModel::DecayWidth_HTotal(){
    return DecayWidth_Hgaga() +DecayWidth_Hee() +DecayWidth_Hmumu() +DecayWidth_Htautau() +DecayWidth_HPiPi() +DecayWidth_HKK() +DecayWidth_HPiPiPiPi() +DecayWidth_Hcc() +DecayWidth_Hss() +DecayWidth_HGluon();
}

// H-> gamma gamma branching fraction
double BSModel::BranchRatio_Hgaga(){
    return DecayWidth_Hgaga()/DecayWidth_HTotal();
}

// H-> e+ e- decay branching fraction
double BSModel::BranchRatio_Hee(){
    return DecayWidth_Hee()/DecayWidth_HTotal();
}

// H-> mu+ mu- decay branching fraction
double BSModel::BranchRatio_Hmumu(){
    return DecayWidth_Hmumu()/DecayWidth_HTotal();
}

// H-> tau+ tau- decay branching fraction
double BSModel::BranchRatio_Htautau(){
    return DecayWidth_Htautau()/DecayWidth_HTotal();
}

// H-> pi pi
double BSModel::BranchRatio_HPiPi(){
    return DecayWidth_HPiPi()/DecayWidth_HTotal();
}

// H-> Kaon Kaon
double BSModel::BranchRatio_HKK(){
    return DecayWidth_HKK()/DecayWidth_HTotal();
}

// H-> pi pi pi pi
double BSModel::BranchRatio_HPiPiPiPi(){
    return DecayWidth_HPiPiPiPi()/DecayWidth_HTotal();
}

// H-> charm charm
double BSModel::BranchRatio_Hcc(){
    return DecayWidth_Hcc()/DecayWidth_HTotal();
}

// H-> strange strange
double BSModel::BranchRatio_Hss(){
    return DecayWidth_Hss()/DecayWidth_HTotal();
}

// H-> gluon gluon branching fraction
double BSModel::BranchRatio_HGluon(){
    return DecayWidth_HGluon()/DecayWidth_HTotal();
}
