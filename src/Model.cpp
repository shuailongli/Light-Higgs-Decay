#include<iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <complex>

# include "JacobiEigenvalue.h"
#include "Model.h"
#include "Parameters.h"
#include "ToolFunctions.h"

using namespace std;


//this function needs to be called to update all parameters whenever changes to input parameters are made
void BSModel::Initialize(){
    
    if(model==0){       //setting all Higgs couplings manuallys
        P11=1, P12=0;
        gHCC=0;
    }
    else if(model==1){  //2HDM Type-I
        xiAuu=xiAcc=xiAtt=1/tanb;
        xiAdd=xiAss=xiAbb=-1/tanb;
        xiAee=xiAmumu=xiAtata=-1/tanb;
        
        xiHuu=xiHcc=xiHtt=cosba-sqrt(1-cosba*cosba)/tanb;
        xiHdd=xiHss=xiHbb=cosba-sqrt(1-cosba*cosba)/tanb;
        xiHee=xiHmumu=xiHtata=cosba-sqrt(1-cosba*cosba)/tanb;
        
        gHCC=-1/vev*(-lambdav*lambdav*(1/tanb-tanb)*sqrt(1-cosba*cosba)+ (2*mHC*mHC-mH*mH+2*lambdav*lambdav)*cosba);
        xiHZZ=cosba, xiHWW=cosba;
        
        P11=1, P12=0;
    }
    else if(model==2){  //2HDM Type-II
        xiAuu=xiAcc=xiAtt=1/tanb;
        xiAdd=xiAss=xiAbb=tanb;
        xiAee=xiAmumu=xiAtata=tanb;
        
        xiHuu=xiHcc=xiHtt=cosba-sqrt(1-cosba*cosba)/tanb;
        xiHdd=xiHss=xiHbb=cosba+sqrt(1-cosba*cosba)*tanb;
        xiHee=xiHmumu=xiHtata=cosba+sqrt(1-cosba*cosba)*tanb;
        
        gHCC=-1/vev*(-lambdav*lambdav*(1/tanb-tanb)*sqrt(1-cosba*cosba)+ (2*mHC*mHC-mH*mH+2*lambdav*lambdav)*cosba);
        xiHZZ=cosba, xiHWW=cosba;
        
        P11=1, P12=0;
    }
    else if(model==3){  //2HDM Type-L
        xiAuu=xiAcc=xiAtt=1/tanb;
        xiAdd=xiAss=xiAbb=-1/tanb;
        xiAee=xiAmumu=xiAtata=tanb;
        
        xiHuu=xiHcc=xiHtt=cosba-sqrt(1-cosba*cosba)/tanb;
        xiHdd=xiHss=xiHbb=cosba-sqrt(1-cosba*cosba)/tanb;
        xiHee=xiHmumu=xiHtata=cosba+sqrt(1-cosba*cosba)*tanb;
        
        gHCC=-1/vev*(-lambdav*lambdav*(1/tanb-tanb)*sqrt(1-cosba*cosba)+ (2*mHC*mHC-mH*mH+2*lambdav*lambdav)*cosba);
        xiHZZ=cosba, xiHWW=cosba;
        
        P11=1, P12=0;
    }
    else if(model==4){  //2HDM Type-F
        xiAuu=xiAcc=xiAtt=1/tanb;
        xiAdd=xiAss=xiAbb=tanb;
        xiAee=xiAmumu=xiAtata=-1/tanb;
        
        xiHuu=xiHcc=xiHtt=cosba-sqrt(1-cosba*cosba)/tanb;
        xiHdd=xiHss=xiHbb=cosba+sqrt(1-cosba*cosba)*tanb;
        xiHee=xiHmumu=xiHtata=cosba-sqrt(1-cosba*cosba)/tanb;
        
        gHCC=-1/vev*(-lambdav*lambdav*(1/tanb-tanb)*sqrt(1-cosba*cosba)+ (2*mHC*mHC-mH*mH+2*lambdav*lambdav)*cosba);
        xiHZZ=cosba, xiHWW=cosba;
        
        P11=1, P12=0;
    }
    else if(model==5){  //MSSM
        xiAuu=xiAcc=xiAtt=1/tanb;
        xiAdd=xiAss=xiAbb=tanb;
        xiAee=xiAmumu=xiAtata=tanb;
        
        xiHuu=xiHcc=xiHtt=cosba-sqrt(1-cosba*cosba)/tanb;
        xiHdd=xiHss=xiHbb=cosba+sqrt(1-cosba*cosba)*tanb;
        xiHee=xiHmumu=xiHtata=cosba+sqrt(1-cosba*cosba)*tanb;
        
        xiHZZ=cosba, xiHWW=cosba;
        
        P11=1, P12=0;
        
    }
    else if(model==6){  //NMSSM
        xiAuu=xiAcc=xiAtt=1/tanb;
        xiAdd=xiAss=xiAbb=tanb;
        xiAee=xiAmumu=xiAtata=tanb;
        
        xiHuu=xiHcc=xiHtt=cosba-sqrt(1-cosba*cosba)/tanb;
        xiHdd=xiHss=xiHbb=cosba+sqrt(1-cosba*cosba)*tanb;
        xiHee=xiHmumu=xiHtata=cosba+sqrt(1-cosba*cosba)*tanb;
        
        xiHZZ=cosba, xiHWW=cosba;
        
        P12=-sqrt(1-P11*P11);
    }
    
    //Higgs effective couplings
    //see eq. (3) in 1612.06538
    //Contribution to Cgamma from quarks and leptons
    Cgamma=-P11/(2*vev)*(3*(2/3.)*(2/3.)*(
            xiAtt*loop_func(4*mt*mt/(mA*mA))            //top
            +xiAcc*loop_func(4*mc*mc/(mA*mA))           //charm
            +xiAuu*loop_func(4*mu*mu/(mA*mA)))          //up
            +3*(-1/3.)*(-1/3.)
            *(xiAbb*loop_func(4*mb*mb/(mA*mA))          //bottom
            +xiAss*loop_func(4*ms*ms/(mA*mA))           //strange
            +xiAdd*loop_func(4*md*md/(mA*mA)))          //down
            +(xiAtata*loop_func(4*mtau*mtau/(mA*mA))    //tau
            +xiAmumu*loop_func(4*mmu*mmu/(mA*mA))       //mu
            +xiAee*loop_func(4*me*me/(mA*mA))));        //e
    
    //Contribution to Cgamma from chargninos
    if(model==5 || model==6){
        
        //Diagonalize the chargino mass matrix for MSSM/NMSSM
        double CharginoMix[2][2]={{WinoMass,mW*sin(atan(tanb))*sqrt(2.)}, {mW*cos(atan(tanb))*sqrt(2.),mueff}};
        
        thetaU=-0.5*atan(2*(CharginoMix[0][0]*CharginoMix[1][0] +CharginoMix[0][1]*CharginoMix[1][1])/(CharginoMix[1][0]*CharginoMix[1][0] +CharginoMix[1][1]*CharginoMix[1][1]-CharginoMix[0][0]*CharginoMix[0][0] -CharginoMix[0][1]*CharginoMix[0][1]));
        thetaV=-0.5*atan(2*(CharginoMix[1][0]*CharginoMix[1][1] +CharginoMix[0][0]*CharginoMix[0][1])/(CharginoMix[0][1]*CharginoMix[0][1] +CharginoMix[1][1]*CharginoMix[1][1]-CharginoMix[0][0]*CharginoMix[0][0] -CharginoMix[1][0]*CharginoMix[1][0]));
        
        //rotation matrices by thetaU and thetaV
        double Umatrix[2][2]={{cos(thetaU),sin(thetaU)},{-sin(thetaU),cos(thetaU)}};
        double Vmatrix[2][2]={{cos(thetaV),sin(thetaV)},{-sin(thetaV),cos(thetaV)}};
        double DiagCharginoMix[2][2]={{0,0},{0,0}};
        for(int i=0;i<2;i++){
            for(int j=0;j<2;j++){
                DiagCharginoMix[i][j]=0;
                for(int k=0;k<2;k++){
                    for(int l=0;l<2;l++){
                        DiagCharginoMix[i][j]+=Umatrix[i][k]*CharginoMix[k][l]*Vmatrix[j][l];
                    }
                }
            }
        }
        
        //in case that mass matrix is not diagonalizd but anti-diagonalized
        //add pi/2 to thetaV to implement diagonalization correctly
        if(abs(DiagCharginoMix[0][0])<abs(DiagCharginoMix[0][1]) && abs(DiagCharginoMix[0][0])<abs(DiagCharginoMix[1][0])){
            thetaV+=PI/2.;
            Vmatrix[0][0]=cos(thetaV), Vmatrix[0][1]=sin(thetaV);
            Vmatrix[1][0]=-sin(thetaV), Vmatrix[1][1]=cos(thetaV);
            
            for(int i=0;i<2;i++){
                for(int j=0;j<2;j++){
                    DiagCharginoMix[i][j]=0;
                    for(int k=0;k<2;k++){
                        for(int l=0;l<2;l++){
                            DiagCharginoMix[i][j]+=Umatrix[i][k]*CharginoMix[k][l]*Vmatrix[j][l];
                        }
                    }
                }
            }
        }
        
        //read chargino masses
        CharginoMass[0]=abs(DiagCharginoMix[0][0]);
        CharginoMass[1]=abs(DiagCharginoMix[1][1]);
        
        //An minus sign is added based on what's given in (3) from 1612.06538 to make the result match.
        for(int i=0;i<2;i++){
            Cgamma+=+Qchargino*Qchargino/(2*sqrt(2.))*(1/CharginoMass[i]*(lambdaSHH*P12*Umatrix[i][1] *Vmatrix[i][1]-2*mW/vev*P11*(cos(atan(tanb))*Umatrix[i][0]*Vmatrix[i][1] +sin(atan(tanb))*Umatrix[i][1]*Vmatrix[i][0])) *loop_func(4*CharginoMass[i]*CharginoMass[i]/(mA*mA)));
        }
    }
    
    
    //Contribution to Cgluon from quarks and leptons
    Cgluon=-P11/(4*vev)*((
            xiAtt*loop_func(4*mt*mt/(mA*mA))        //top
            +xiAcc*loop_func(4*mc*mc/(mA*mA))       //charm
            +xiAuu*loop_func(4*mu*mu/(mA*mA)))      //up
            +(xiAbb*loop_func(4*mb*mb/(mA*mA))      //bottom
            +xiAss*loop_func(4*ms*ms/(mA*mA))       //strange
            +xiAdd*loop_func(4*md*md/(mA*mA))));    //down
    
    
    //effective Yukawa couplings in spectator model
    YuAu=sqrt(2.)*P11/(sqrt(3.)*vev*fpi)*Bmudf*xiAuu;
    YuAd=sqrt(2.)*P11/(sqrt(3.)*vev*fpi)*Bmddf*xiAdd;
    YuAs=sqrt(2.)*P11/(sqrt(3.)*vev*fpi)*Bmsdf*xiAss;
    
    
    //A-meson mixing matrix
    double NCgluon=sqrt(real(Cgluon*conj(Cgluon)));
    double dm32=-fpi*P11/(2*vev)*(Bmudf*xiAuu-Bmddf*xiAdd);
    double dm82=-fpi*P11/(2*sqrt(3.)*vev)*(Bmudf*xiAuu+Bmddf*xiAdd-2*Bmsdf*xiAss);
    double dm92=-fpi*P11/(sqrt(6.)*vev)*(Bmudf*xiAuu+Bmddf*xiAdd+Bmsdf*xiAss) -sqrt(3/2.)*fpi*NCgluon*invCf2;
    double barmA2=mA*mA+NCgluon*NCgluon*fpi*fpi*invCf2;
    double MassMix[4][4]={{mpi*mpi,0,0,dm32},
        {0,mpi8*mpi8,Delta,dm82},
        {0,Delta,mpi9*mpi9,dm92},
        {dm32,dm82,dm92,barmA2}};
    
    
    double varMatrix[16];
    double eigenvec[16], eigenval[4];
    int it_num, rot_num;
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++) varMatrix[j+i*4]=MassMix[i][j];
    }
    jacobi_eigenvalue(4, varMatrix, 100, eigenvec, eigenval, it_num, rot_num);
    
    double eigenA[4];
    for(int i=0;i<16;i++){
        eigenA[i%4]=eigenvec[i];
        if((i+1)%4==0){
            if(eigenA[3]>eigenA[2] && eigenA[3]>eigenA[1] && eigenA[3]>eigenA[0]) break;
        }
    }
    
    double mixA8=eigenA[1], mixA9=eigenA[2];
    mixAA=eigenA[3];
    mixA3=eigenA[0];
    mixAeta=coseta*mixA8-sineta*mixA9;
    mixAetap=sineta*mixA8+coseta*mixA9;
    
}


//initialize xiAtt, xiAbb, xiAcc, xiAss, xiAuu, xiAdd, xiAtata, xiAmumu, xiAee
//can only be invoked for model==0
void BSModel::set_AffCoupling(double xit,double xib,double xic,double xis,double xiu,double xid,double xitau,double ximu,double xie){
    if(model!=0){
        cout<<"*Error: set_AffCoupling() is only used for model=0."<<endl;
        cout<<"        Please set model to 0."<<endl;
        exit(0);
    }
    else{
        xiAtt=xit, xiAbb=xib, xiAcc=xic;
        xiAss=xis, xiAuu=xiu, xiAdd=xid;
        xiAtata=xitau, xiAmumu=ximu, xiAee=xie;
    }
}

//initialize xiHtt, xiHbb, xiHcc, xiHss, xiHuu, xiHdd, xiHtata, xiHmumu, xiHee
//can only be invoked for model==0
void BSModel::set_HffCoupling(double xit,double xib,double xic,double xis,double xiu,double xid,double xitau,double ximu,double xie){
    if(model!=0){
        cout<<"*Error: set_HffCoupling() is only used for model=0."<<endl;
        cout<<"        Please set model=0 before use it."<<endl;
        exit(0);
    }
    else{
        xiHtt=xit, xiHbb=xib, xiHcc=xic;
        xiHss=xis, xiHuu=xiu, xiHdd=xid;
        xiHtata=xitau, xiHmumu=ximu, xiHee=xie;
    }
}

//initialize xiHZZ, xiHWW
//can only be invoked for model==0
void BSModel::set_HVVCoupling(double xiz,double xiw){
    if(model!=0){
        cout<<"*Error: set_HVVCoupling() is only used for model=0."<<endl;
        cout<<"        Please set model=0 before use it."<<endl;
        exit(0);
    }
    else{
        xiHZZ=xiz, xiHWW=xiw;
    }
}

//initialize gHCC
//can only be invoked for model==0
void BSModel::set_gHCC(double HCC){
    if(model!=0){
        cout<<"*Error: set_gHCC() is only used for model=0."<<endl;
        cout<<"        Please set model=0 before use it."<<endl;
        exit(0);
    }
    else{
        gHCC=HCC;
    }
}


// A-> gamma gamma decay width
double BSModel::DecayWidth_Agaga(){
    //effective coupling of pi0 gamma gamma: -10.75 GeV^-1
    //effective coupling of eta gamma gamma: -10.8 GeV^-1
    //effective coupling of eta' gamma gamma: -13.6 GeV^-1
    //see eq. (19) in 1612.06538
    complex<double> sumA=-mixAA*Cgamma+mixA3*(-10.75) +mixAeta*(-10.8)+mixAetap*(-13.6);
    return real((1/(137.*137))*mA*mA*mA/(64*PI*PI *PI)*sumA*conj(sumA));
}

// A-> e+ e- decay width
double BSModel::DecayWidth_Aee(){
    //pure Higgs state decay dominates, no mixing included
    double YAee=me/vev*P11*xiAee;
    return (1-4*me*me/(mA*mA)>0)? YAee*YAee*mA/(8*PI)*sqrt(1-4*me*me/(mA*mA)):0;
}

// A-> mu+ mu- decay width
double BSModel::DecayWidth_Amumu(){
    //pure Higgs state decay dominates, no mixing included
    double YAmm=mmu/vev*P11*xiAmumu;
    return (1-4*mmu*mmu/(mA*mA)>0)? YAmm*YAmm*mA/(8*PI)*sqrt(1-4*mmu*mmu/(mA*mA)):0;
}


// A-> tau+ tau- decay width
double BSModel::DecayWidth_Atautau(){
    //pure Higgs state decay dominates, no mixing included
    double YAtautau=mtau/vev*P11*xiAtata;
    return (1-4*mtau*mtau/(mA*mA)>0)? YAtautau*YAtautau*mA/(8*PI)*sqrt(1-4*mtau*mtau/(mA*mA)):0;
}


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
double BSModel::cupAPiPiPi(int pi1, int pi2, int pi3){
    double product[3][3];
    PermuProduct(Meson[pi1-1], Meson[pi2-1], Meson[pi3-1], product);
    return -sqrt(2.)*P11/(6.*vev*fpi)*(product[0][0]*Bmudf*xiAuu + product[1][1]*Bmddf*xiAdd + product[2][2]*Bmsdf*xiAss);
}

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
double BSModel::cupPiPiPiPi(int pi1, int pi2, int pi3, int pi4){
    //sort pi1, pi2, pi3, pi4 according to ascending order
    int piID[4]={pi1, pi2, pi3, pi4};
    for(int i=0;i<4;i++){
        for(int j=i;j<4;j++){
            if(piID[i]>piID[j]){
                int temp=piID[i];
                piID[i]=piID[j];
                piID[j]=temp;
            }
        }
    }
    
    double tmp_delta=2*(mKpm*mKpm-mK0*mK0-0.1395*0.1395+0.135*0.135); //tmp_delta=Bmudf-Bmddf
    
    if(piID[3]==8 && piID[2]==3 && piID[1]==3 && piID[0]==3) return 0.721729;
    else if(piID[3]==8 && piID[2]==3 && piID[1]==2 && piID[0]==1) return 0.246559;
    else if(piID[3]==9 && piID[2]==3 && piID[1]==3 && piID[0]==3) return 0.299175;
    else if(piID[3]==9 && piID[2]==3 && piID[1]==2 && piID[0]==1) return 0.145544;
    else if(piID[3]==9 && piID[2]==8 && piID[1]==3 && piID[0]==3) return 6.40;
    else if(piID[3]==9 && piID[2]==8 && piID[1]==2 && piID[0]==1) return 6.1967;
    else if(piID[3]==9 && piID[2]==9 && piID[1]==3 && piID[0]==3) return 5.70933;
    else if(piID[3]==9 && piID[2]==9 && piID[1]==2 && piID[0]==1) return 5.52797;
    else if(piID[3]==8 && piID[2]==8 && piID[1]==8 && piID[0]==3){
        return 1/(6*sqrt(3.)*fpi*fpi)*tmp_delta* (coseta-sqrt(2.)*sineta) *(coseta-sqrt(2.)*sineta) *(coseta-sqrt(2.)*sineta);
    }
    else if(piID[3]==9 && piID[2]==8 && piID[1]==8 && piID[0]==3){
        return 1/(6*sqrt(3.)*fpi*fpi)*tmp_delta* (coseta-sqrt(2.)*sineta) *(coseta-sqrt(2.)*sineta) *(sineta+sqrt(2.)*coseta);
    }
    else if(piID[3]==9 && piID[2]==9 && piID[1]==8 && piID[0]==3){
        return 1/(6*sqrt(3.)*fpi*fpi)*tmp_delta* (coseta-sqrt(2.)*sineta) *(sineta+sqrt(2.)*coseta) *(sineta+sqrt(2.)*coseta);
    }
    else if(piID[3]==9 && piID[2]==9 && piID[1]==9 && piID[0]==3){
        return 1/(6*sqrt(3.)*fpi*fpi)*tmp_delta* (sineta+sqrt(2.)*coseta) *(sineta+sqrt(2.)*coseta) *(sineta+sqrt(2.)*coseta);
    }
    else{
        double product[3][3];
        PermuProduct(Meson[pi1-1], Meson[pi2-1], Meson[pi3-1], Meson[pi4-1], product);
        return 1/(12.*fpi*fpi)*(product[0][0]*Bmudf + product[1][1]*Bmddf + product[2][2]*Bmsdf);
    }
}



// A-> pi pi pi
double BSModel::DecayWidth_APiPiPi(){
    double cupMixingP0P0P0; // A-> pi0 pi0 pi0
    double cupMixingP0PpPm; // A-> pi0 pi+ pi-
    if(mA>1.3){
        double AA_PiPiPi=5*3./144*(YuAu*YuAu+YuAd*YuAd);
        double cup_AP0P0P0=sqrt(18/5.*AA_PiPiPi)*SIGN(cupAPiPiPi(3,3,3));
        double cup_AP0PpPm=cup_AP0P0P0/3.;
        cupMixingP0P0P0=cup_AP0P0P0*mixAA +cupPiPiPiPi(3,3,3,3)*mixA3 +cupPiPiPiPi(8,3,3,3)*mixAeta +cupPiPiPiPi(9,3,3,3)*mixAetap;
        cupMixingP0PpPm=cup_AP0PpPm*mixAA +cupPiPiPiPi(3,3,2,1)*mixA3 +cupPiPiPiPi(8,3,2,1)*mixAeta +cupPiPiPiPi(9,3,2,1)*mixAetap;
    }
    else{
        
        cupMixingP0P0P0=cupAPiPiPi(3,3,3)*mixAA +cupPiPiPiPi(3,3,3,3)*mixA3 +cupPiPiPiPi(8,3,3,3)*mixAeta +cupPiPiPiPi(9,3,3,3)*mixAetap;
        cupMixingP0PpPm=cupAPiPiPi(3,2,1)*mixAA +cupPiPiPiPi(3,3,2,1)*mixA3 +cupPiPiPiPi(8,3,2,1)*mixAeta +cupPiPiPiPi(9,3,2,1)*mixAetap;
    }
    return 1/(256.*PI*PI*PI*mA) *(cupMixingP0P0P0*cupMixingP0P0P0/6. +cupMixingP0PpPm*cupMixingP0PpPm)*APPP(mA, mpi, mpi, mpi);
}


// A-> eta pi pi
double BSModel::DecayWidth_AEtaPiPi(){
    double cupMixingEtaP0P0;    // A-> eta pi0 pi0
    double cupMixingEtaPpPm;    // A-> eta pi+ pi-
    if(mA>1.3){
        double AA_EtaPiPi=3./16*(YuAu*YuAu+YuAd*YuAd) *(coseta-sqrt(2.)*sineta) *(coseta-sqrt(2.)*sineta);
        double cup_AEtaP0P0=sqrt(2/3.*AA_EtaPiPi)*SIGN(cupAPiPiPi(8,3,3));
        double cup_AEtaPpPm=cup_AEtaP0P0;
        cupMixingEtaP0P0=cup_AEtaP0P0*mixAA +cupPiPiPiPi(8,3,3,3)*mixA3 +cupPiPiPiPi(8,8,3,3)*mixAeta +cupPiPiPiPi(9,8,3,3)*mixAetap;
        cupMixingEtaPpPm=cup_AEtaPpPm*mixAA +cupPiPiPiPi(8,3,2,1)*mixA3 +cupPiPiPiPi(8,8,2,1)*mixAeta +cupPiPiPiPi(9,8,2,1)*mixAetap;
    }
    else{
        cupMixingEtaP0P0=cupAPiPiPi(8,3,3)*mixAA +cupPiPiPiPi(8,3,3,3)*mixA3 +cupPiPiPiPi(8,8,3,3)*mixAeta +cupPiPiPiPi(9,8,3,3)*mixAetap;
        cupMixingEtaPpPm=cupAPiPiPi(8,2,1)*mixAA +cupPiPiPiPi(8,3,2,1)*mixA3 +cupPiPiPiPi(8,8,2,1)*mixAeta +cupPiPiPiPi(9,8,2,1)*mixAetap;
    }
    return 1/(256.*PI*PI*PI*mA) *(cupMixingEtaP0P0*cupMixingEtaP0P0/2. +cupMixingEtaPpPm*cupMixingEtaPpPm)*APPP(mA, meta, mpi, mpi);
}


// A-> eta' pi pi
double BSModel::DecayWidth_AEtapPiPi(){
    double cupMixingEtapP0P0;   // A-> eta' pi0 pi0
    double cupMixingEtapPpPm;   // A-> eta' pi+ pi-
    if(mA>1.3){
        double AA_EtapPiPi=3./16*(YuAu*YuAu+YuAd*YuAd) *(sineta+sqrt(2.)*coseta) *(sineta+sqrt(2.)*coseta);
        double cup_AEtapP0P0=sqrt(2/3.*AA_EtapPiPi)*SIGN(cupAPiPiPi(9,3,3));
        double cup_AEtapPpPm=cup_AEtapP0P0;
        cupMixingEtapP0P0=cup_AEtapP0P0*mixAA +cupPiPiPiPi(9,3,3,3)*mixA3 +cupPiPiPiPi(9,8,3,3)*mixAetap +cupPiPiPiPi(9,9,3,3)*mixAeta;
        cupMixingEtapPpPm=cup_AEtapPpPm*mixAA +cupPiPiPiPi(9,3,2,1)*mixA3 +cupPiPiPiPi(9,8,2,1)*mixAetap +cupPiPiPiPi(9,9,2,1)*mixAeta;
    }
    else{
        //interchanged mixAetap and mixAeta to match arXiv: 1612.06538
        cupMixingEtapP0P0=cupAPiPiPi(9,3,3)*mixAA +cupPiPiPiPi(9,3,3,3)*mixA3 +cupPiPiPiPi(9,8,3,3)*mixAetap +cupPiPiPiPi(9,9,3,3)*mixAeta;
        cupMixingEtapPpPm=cupAPiPiPi(9,2,1)*mixAA +cupPiPiPiPi(9,3,2,1)*mixA3 +cupPiPiPiPi(9,8,2,1)*mixAetap +cupPiPiPiPi(9,9,2,1)*mixAeta;
    }
    
    
    return 1/(256.*PI*PI*PI*mA) *(cupMixingEtapP0P0*cupMixingEtapP0P0/2.+cupMixingEtapPpPm*cupMixingEtapPpPm) *APPP(mA, metap, mpi, mpi);
    
}

// A-> eta eta pi
double BSModel::DecayWidth_AEtaEtaPi(){
    double cupMixingP0EtaEta;
    if(mA>1.3){
        double AA_PiEtaEta=3./144*(YuAu*YuAu+YuAd*YuAd) *pow(coseta-sqrt(2.)*sineta,4);
        double cup_AEtaEtaP0=sqrt(AA_PiEtaEta*2.)*SIGN(cupAPiPiPi(3,8,8));
        cupMixingP0EtaEta=cupAPiPiPi(3,8,8)*mixAA +cupPiPiPiPi(8,8,3,3)*mixA3 +cupPiPiPiPi(8,8,8,3)*mixAeta +cupPiPiPiPi(9,8,8,3)*mixAetap;
    }
    else{
        cupMixingP0EtaEta=cupAPiPiPi(3,8,8)*mixAA +cupPiPiPiPi(8,8,3,3)*mixA3 +cupPiPiPiPi(8,8,8,3)*mixAeta +cupPiPiPiPi(9,8,8,3)*mixAetap;
    }
    
    return 1/(256.*PI*PI*PI*mA) *(cupMixingP0EtaEta*cupMixingP0EtaEta/2. ) *APPP(mA, meta, meta, mpi);
}


// A-> K K pi
double BSModel::DecayWidth_AKKPi(){
    double cupMixingP0KpKm;     // A-> pi0 K+ k-
    double cupMixingP0K0BK0;    // A-> pi0 K0 BK0
    double cupMixingPpKmK0;     // A-> pi+ K- K0
    double cupMixingPmKpBK0;    // A-> pi- K+ BK0
    if(mA>1.3){
        double AA_PiKK=3./36*(4*YuAu*YuAu+4*YuAd*YuAd+3*YuAs*YuAs);
        double cup_AP0KpKm=sqrt(3/36.*(2*YuAu*YuAu+1/2.*YuAs*YuAs))*SIGN(cupAPiPiPi(3,4,5));
        double cup_AP0K0BK0=sqrt(3/36.*(2*YuAd*YuAd+1/2.*YuAs*YuAs))*SIGN(cupAPiPiPi(3,6,7));
        double cup_APpKmK0=sqrt(3/36.*(YuAu*YuAu+YuAd*YuAd+YuAs*YuAs))*SIGN(cupAPiPiPi(1,5,6));
        double cup_APmKpBK0=sqrt(3/36.*(YuAu*YuAu+YuAd*YuAd+YuAs*YuAs))*SIGN(cupAPiPiPi(2,4,7));
        cupMixingP0KpKm=cup_AP0KpKm*mixAA +cupPiPiPiPi(3,3,4,5)*mixA3 +cupPiPiPiPi(8,3,4,5)*mixAeta +cupPiPiPiPi(9,3,4,5)*mixAetap;
        cupMixingP0K0BK0=cup_AP0K0BK0*mixAA +cupPiPiPiPi(3,3,6,7)*mixA3 +cupPiPiPiPi(8,3,6,7)*mixAeta +cupPiPiPiPi(9,3,6,7)*mixAetap;
        cupMixingPpKmK0=cup_APpKmK0*mixAA +cupPiPiPiPi(3,1,5,6)*mixA3 +cupPiPiPiPi(8,1,5,6)*mixAeta +cupPiPiPiPi(9,1,5,6)*mixAetap;
        cupMixingPmKpBK0=cup_APmKpBK0*mixAA +cupPiPiPiPi(3,2,4,7)*mixA3 +cupPiPiPiPi(8,2,4,7)*mixAeta +cupPiPiPiPi(9,2,4,7)*mixAetap;
    }
    else{
        cupMixingP0KpKm=cupAPiPiPi(3,4,5)*mixAA +cupPiPiPiPi(3,3,4,5)*mixA3 +cupPiPiPiPi(8,3,4,5)*mixAeta +cupPiPiPiPi(9,3,4,5)*mixAetap;
        cupMixingP0K0BK0=cupAPiPiPi(3,6,7)*mixAA +cupPiPiPiPi(3,3,6,7)*mixA3 +cupPiPiPiPi(8,3,6,7)*mixAeta +cupPiPiPiPi(9,3,6,7)*mixAetap;
        cupMixingPpKmK0=cupAPiPiPi(1,5,6)*mixAA +cupPiPiPiPi(3,1,5,6)*mixA3 +cupPiPiPiPi(8,1,5,6)*mixAeta +cupPiPiPiPi(9,1,5,6)*mixAetap;
        cupMixingPmKpBK0=cupAPiPiPi(2,4,7)*mixAA +cupPiPiPiPi(3,2,4,7)*mixA3 +cupPiPiPiPi(8,2,4,7)*mixAeta +cupPiPiPiPi(9,2,4,7)*mixAetap;
    }
    
    return 1/(256.*PI*PI*PI*mA) *(cupMixingP0KpKm*cupMixingP0KpKm  +cupMixingP0K0BK0*cupMixingP0K0BK0 +cupMixingPpKmK0*cupMixingPpKmK0 +cupMixingPmKpBK0*cupMixingPmKpBK0) *APPP(mA, mpi, mKpm, mKpm);
}


// A-> gamma pi pi
double BSModel::DecayWidth_AGammaPiPi(){
    double a=4*mpi*mpi, b=mA*mA;
    if(a>=b) return 0;
    else{
        double dx=(b-a)/10000;
        double sum=0;
        for(double x = a; x < b; x += dx)
        {
            
            double func_x=Gamma0(mA,x)*real((mixAeta *BEta(mA,x).first + mixAetap *BEta(mA,x).second)*conj(mixAeta *BEta(mA,x).first + mixAetap *BEta(mA,x).second));
            double func_x_dxd2=Gamma0(mA,x + dx/2)*real((mixAeta *BEta(mA,x + dx/2).first + mixAetap *BEta(mA,x + dx/2).second)*conj(mixAeta *BEta(mA,x + dx/2).first + mixAetap *BEta(mA,x + dx/2).second));
            double func_x_dx=Gamma0(mA,x + dx)*real((mixAeta *BEta(mA,x + dx).first + mixAetap *BEta(mA,x + dx).second)*conj(mixAeta *BEta(mA,x + dx).first + mixAetap *BEta(mA,x + dx).second));
            sum += func_x + 4*func_x_dxd2 + func_x_dx;
        }
        
        sum *= dx/6;
        
        return sum;
    }
}



// A-> eta eta' pi
double BSModel::DecayWidth_AEtaEtapPi(){
    double cupMixingEtaEtapP0;
    if(mA<1.3){
        double AA_PiEtaEtap=3./72*(YuAu*YuAu+YuAd*YuAd) *pow(coseta-sqrt(2.)*sineta,2) *pow(sineta+sqrt(2.)*coseta,2);
        double cup_AEtaEtapP0=sqrt(AA_PiEtaEtap)*SIGN(cupAPiPiPi(3,8,9));
        cupMixingEtaEtapP0=cup_AEtaEtapP0*mixAA +cupPiPiPiPi(3,3,8,9)*mixA3 +cupPiPiPiPi(8,3,8,9)*mixAeta +cupPiPiPiPi(9,3,8,9)*mixAetap;
    }
    else{
        cupMixingEtaEtapP0=cupAPiPiPi(3,8,9)*mixAA +cupPiPiPiPi(3,3,8,9)*mixA3 +cupPiPiPiPi(8,3,8,9)*mixAeta +cupPiPiPiPi(9,3,8,9)*mixAetap;
    }
    return 1/(256.*PI*PI*PI*mA) *(cupMixingEtaEtapP0*cupMixingEtaEtapP0  ) *APPP(mA, mpi, meta, metap);
}


// A-> eta' eta' pi
double BSModel::DecayWidth_AEtapEtapPi(){
    double cupEtapEtapP0;
    if(mA>1.3){
        double AA_PiEtapEtap=3./144*(YuAu*YuAu+YuAd*YuAd) *pow(sineta+sqrt(2.)*coseta,4);
        double cupAEtapEtapP0=sqrt(AA_PiEtapEtap*2.)*SIGN(cupAPiPiPi(3,9,9));
        cupEtapEtapP0=cupAEtapEtapP0*mixAA +cupPiPiPiPi(3,3,9,9)*mixA3 +cupPiPiPiPi(8,3,9,9)*mixAeta +cupPiPiPiPi(9,3,9,9)*mixAetap;
    }
    else{
        cupEtapEtapP0=cupAPiPiPi(3,9,9)*mixAA +cupPiPiPiPi(3,3,9,9)*mixA3 +cupPiPiPiPi(8,3,9,9)*mixAeta +cupPiPiPiPi(9,3,9,9)*mixAetap;
    }
    return 1/(256.*PI*PI*PI*mA) *(cupEtapEtapP0*cupEtapEtapP0/2.  ) *APPP(mA, mpi, metap, metap);
}


// A-> eta eta eta
double BSModel::DecayWidth_AEtaEtaEta(){
    double cupEtaEtaEta;
    if(mA>1.3){
        double AA_EtaEtaEta=3./1296*(YuAu*YuAu+YuAd*YuAd+64*YuAs*YuAs);
        double cup_AEtaEtaEta=sqrt(6.*AA_EtaEtaEta)*SIGN(cupAPiPiPi(8,8,8));
        cupEtaEtaEta=cup_AEtaEtaEta*mixAA +cupPiPiPiPi(3,8,8,8)*mixA3 +cupPiPiPiPi(8,8,8,8)*mixAeta +cupPiPiPiPi(9,8,8,8)*mixAetap;
    }
    else{
        cupEtaEtaEta=cupAPiPiPi(8,8,8)*mixAA +cupPiPiPiPi(3,8,8,8)*mixA3 +cupPiPiPiPi(8,8,8,8)*mixAeta +cupPiPiPiPi(9,8,8,8)*mixAetap;
    }
    return 1/(256.*PI*PI*PI*mA) *(cupEtaEtaEta*cupEtaEtaEta/6.  ) *APPP(mA, meta, meta, meta);
}


// A-> eta eta eta'
double BSModel::DecayWidth_AEtaEtaEtap(){
    double cupEtaEtaEtap;
    if(mA>1.3){
        double AA_EtaEtaEtap=3./216*(YuAu*YuAu+YuAd*YuAd+16*YuAs*YuAs);
        double cup_AEtaEtaEtap=sqrt(2.*AA_EtaEtaEtap)*SIGN(cupAPiPiPi(8,8,9));
        cupEtaEtaEtap=cup_AEtaEtaEtap*mixAA +cupPiPiPiPi(3,8,8,9)*mixA3 +cupPiPiPiPi(8,8,8,9)*mixAeta +cupPiPiPiPi(9,8,8,9)*mixAetap;
    }
    else{
        cupEtaEtaEtap=cupAPiPiPi(8,8,9)*mixAA +cupPiPiPiPi(3,8,8,9)*mixA3 +cupPiPiPiPi(8,8,8,9)*mixAeta +cupPiPiPiPi(9,8,8,9)*mixAetap;
    }
    return 1/(256.*PI*PI*PI*mA) *(cupEtaEtaEtap*cupEtaEtaEtap/2.  ) *APPP(mA, meta, meta, metap);
}


// A-> eta eta' eta'
double BSModel::DecayWidth_AEtaEtapEtap(){
    double cupEtaEtapEtap;
    if(mA>1.3){
        double AA_EtaEtapEtap=3./108*(YuAu*YuAu+YuAd*YuAd+4*YuAs*YuAs);
        double cup_AEtaEtapEtap=sqrt(2.*AA_EtaEtapEtap)*SIGN(cupAPiPiPi(8,9,9));
        cupEtaEtapEtap=cup_AEtaEtapEtap*mixAA +cupPiPiPiPi(3,8,9,9)*mixA3 +cupPiPiPiPi(8,8,9,9)*mixAeta +cupPiPiPiPi(9,8,9,9)*mixAetap;
    }
    else{
        cupEtaEtapEtap=cupAPiPiPi(8,9,9)*mixAA +cupPiPiPiPi(3,8,9,9)*mixA3 +cupPiPiPiPi(8,8,9,9)*mixAeta +cupPiPiPiPi(9,8,9,9)*mixAetap;
    }
    return 1/(256.*PI*PI*PI*mA) *(cupEtaEtapEtap*cupEtaEtapEtap/2.  ) *APPP(mA, meta, metap, metap);
}


// A-> eta' eta' eta'
double BSModel::DecayWidth_AEtapEtapEtap(){
    double cupEtapEtapEtap;
    if(mA>1.3){
        double AA_EtapEtapEtap=3./162*(YuAu*YuAu+YuAd*YuAd+YuAs*YuAs);
        double cup_EtapEtapEtap=sqrt(6.*AA_EtapEtapEtap)*SIGN(cupAPiPiPi(9,9,9));
        cupEtapEtapEtap=cup_EtapEtapEtap*mixAA +cupPiPiPiPi(3,9,9,9)*mixA3 +cupPiPiPiPi(8,9,9,9)*mixAeta +cupPiPiPiPi(9,9,9,9)*mixAetap;
    }
    else{
        cupEtapEtapEtap=cupAPiPiPi(9,9,9)*mixAA +cupPiPiPiPi(3,9,9,9)*mixA3 +cupPiPiPiPi(8,9,9,9)*mixAeta +cupPiPiPiPi(9,9,9,9)*mixAetap;
    }
    return 1/(256.*PI*PI*PI*mA) *(cupEtapEtapEtap*cupEtapEtapEtap/6.  ) *APPP(mA, metap, metap, metap);
}


// A-> eta K K
double BSModel::DecayWidth_AEtaKK(){
    double cupEtaK0BK0; // A-> eta K0 BK0
    double cupEtaKpKm;  // A-> eta K+ K-
    if(mA>1.3){
        double AA_EtaKK=3./12*(YuAs*YuAs);
        double cup_AEtaK0BK0=sqrt(AA_EtaKK/2.)*SIGN(cupAPiPiPi(8,6,7));
        double cup_AEtaKpKm=cup_AEtaK0BK0;
        cupEtaK0BK0=cup_AEtaK0BK0*mixAA +cupPiPiPiPi(3,8,6,7)*mixA3 +cupPiPiPiPi(8,8,6,7)*mixAeta +cupPiPiPiPi(9,8,6,7)*mixAetap;
        cupEtaKpKm=cup_AEtaKpKm*mixAA +cupPiPiPiPi(3,8,4,5)*mixA3 +cupPiPiPiPi(8,8,4,5)*mixAeta +cupPiPiPiPi(9,8,4,5)*mixAetap;
    }
    else{
        cupEtaK0BK0=cupAPiPiPi(8,6,7)*mixAA +cupPiPiPiPi(3,8,6,7)*mixA3 +cupPiPiPiPi(8,8,6,7)*mixAeta +cupPiPiPiPi(9,8,6,7)*mixAetap;
        cupEtaKpKm=cupAPiPiPi(8,4,5)*mixAA +cupPiPiPiPi(3,8,4,5)*mixA3 +cupPiPiPiPi(8,8,4,5)*mixAeta +cupPiPiPiPi(9,8,4,5)*mixAetap;
    }
    return 1/(256.*PI*PI*PI*mA) *(cupEtaK0BK0*cupEtaK0BK0 + cupEtaKpKm*cupEtaKpKm ) *APPP(mA, meta, mK0, mK0);
}


// A-> eta' K K
double BSModel::DecayWidth_AEtapKK(){
    double cupEtapK0BK0;    // A-> eta' K0 BK0
    double cupEtapKpKm;     // A-> eta' K+ K-
    if(mA>1.3){
        double AA_EtapKK=3./12*(YuAu*YuAu+YuAd*YuAd+2*YuAs*YuAs);
        double cup_EtapK0BK0=sqrt(3./12*(YuAd*YuAd+YuAs*YuAs)) *SIGN(cupAPiPiPi(9,6,7));
        double cup_EtapKpKm=sqrt(3./12*(YuAu*YuAu+YuAs*YuAs)) *SIGN(cupAPiPiPi(9,4,5));
        cupEtapK0BK0=cup_EtapK0BK0*mixAA +cupPiPiPiPi(3,9,6,7)*mixA3 +cupPiPiPiPi(8,9,6,7)*mixAeta +cupPiPiPiPi(9,9,6,7)*mixAetap;
        cupEtapKpKm=cup_EtapKpKm*mixAA +cupPiPiPiPi(3,9,4,5)*mixA3 +cupPiPiPiPi(8,9,4,5)*mixAeta +cupPiPiPiPi(9,9,4,5)*mixAetap;
    }
    else{
        cupEtapK0BK0=cupAPiPiPi(9,6,7)*mixAA +cupPiPiPiPi(3,9,6,7)*mixA3 +cupPiPiPiPi(8,9,6,7)*mixAeta +cupPiPiPiPi(9,9,6,7)*mixAetap;
        cupEtapKpKm=cupAPiPiPi(9,4,5)*mixAA +cupPiPiPiPi(3,9,4,5)*mixA3 +cupPiPiPiPi(8,9,4,5)*mixAeta +cupPiPiPiPi(9,9,4,5)*mixAetap;
    }
    return 1/(256.*PI*PI*PI*mA) *(cupEtapK0BK0*cupEtapK0BK0 + cupEtapKpKm*cupEtapKpKm ) *APPP(mA, metap, mK0, mK0);
    }

//partonic hardon decay, mA>3 GeV
// A-> quark quark
double BSModel::DecayWidth_AQuark(){
    double YAuu=mu*P11/vev*xiAuu;
    double YAdd=md*P11/vev*xiAdd;
    double YAss=ms*P11/vev*xiAss;
    double YAcc=mc*P11/vev*xiAcc;
    double YAbb=mb*P11/vev*xiAbb;
    double YAtt=mt*P11/vev*xiAtt;
    
    double Ga_Auu=(1-4*mu*mu/(mA*mA)>0)? 3*YAuu*YAuu*mA/(8*PI)*sqrt(1-4*mu*mu/(mA*mA)):0;
    double Ga_Add=(1-4*md*md/(mA*mA)>0)? 3*YAdd*YAdd*mA/(8*PI)*sqrt(1-4*md*md/(mA*mA)):0;
    double Ga_Ass=(1-4*ms*ms/(mA*mA)>0)? 3*YAss*YAss*mA/(8*PI)*sqrt(1-4*ms*ms/(mA*mA)):0;
    double Ga_Acc=(1-4*mc*mc/(mA*mA)>0)? 3*YAcc*YAcc*mA/(8*PI)*sqrt(1-4*mc*mc/(mA*mA)):0;
    double Ga_Abb=(1-4*mb*mb/(mA*mA)>0)? 3*YAbb*YAbb*mA/(8*PI)*sqrt(1-4*mb*mb/(mA*mA)):0;
    double Ga_Att=(1-4*mt*mt/(mA*mA)>0)? 3*YAtt*YAtt*mA/(8*PI)*sqrt(1-4*mt*mt/(mA*mA)):0;
    
    if(mA<3) return 0;
    else return Ga_Auu+Ga_Add+Ga_Ass+Ga_Acc+Ga_Abb+Ga_Att;
}


// A-> gluon gluon
double BSModel::DecayWidth_AGluon(){
    if(mA<3) return 0;
    else return real((alphasQ(mA)*alphasQ(mA))*mA*mA*mA/(8*PI*PI*PI)*Cgluon*conj(Cgluon));
}


//Total hadronic decay width
double BSModel::DecayWidth_ATotalHadron(){
    if(mA<3){
        return DecayWidth_APiPiPi()
        +DecayWidth_AEtaPiPi()
        +DecayWidth_AEtapPiPi()
        +DecayWidth_AEtaEtaPi()
        +DecayWidth_AKKPi()
        +DecayWidth_AGammaPiPi()
        +DecayWidth_AEtaEtapPi()
        +DecayWidth_AEtapEtapPi()
        +DecayWidth_AEtaEtaEta()
        +DecayWidth_AEtaEtaEtap()
        +DecayWidth_AEtaEtapEtap()
        +DecayWidth_AEtapEtapEtap()
        +DecayWidth_AEtaKK()
        +DecayWidth_AEtapKK();
    }
    else{
        return DecayWidth_AQuark()+DecayWidth_AGluon();
    }
    
}


//Total leptonic decay width
double BSModel::DecayWidth_ATotalLepton(){
    return DecayWidth_Aee()+DecayWidth_Amumu()+DecayWidth_Atautau();
}


//Total decay width
double BSModel::DecayWidth_ATotal(){
    return DecayWidth_ATotalHadron()+DecayWidth_ATotalLepton()+DecayWidth_Agaga();
}



//Decay Branching Fractions
// A-> gamma gamma branching fraction
double BSModel::BranchRatio_Agaga(){
    return DecayWidth_Agaga()/DecayWidth_ATotal();
}


//leptonic decay Branching Fractions
// A-> e+ e- decay branching fraction
double BSModel::BranchRatio_Aee(){
    return DecayWidth_Aee()/DecayWidth_ATotal();
}

// A-> mu+ mu- decay branching fraction
double BSModel::BranchRatio_Amumu(){
    return DecayWidth_Amumu()/DecayWidth_ATotal();
}

// A-> tau+ tau- decay branching fraction
double BSModel::BranchRatio_Atautau(){
    return DecayWidth_Atautau()/DecayWidth_ATotal();
}

//Total hadronic decay Branching Fractions
double BSModel::BranchRatio_ATotalHadron(){
    return DecayWidth_ATotalHadron()/DecayWidth_ATotal();
}

//partonic hardon decay, mA>3 GeV
// A-> quark quark branching fraction
double BSModel::BranchRatio_AQuark(){
    return DecayWidth_AQuark()/DecayWidth_ATotal();
}

// A-> gluon gluon branching fraction
double BSModel::BranchRatio_AGluon(){
    return DecayWidth_AGluon()/DecayWidth_ATotal();
}
