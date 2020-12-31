#include<iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <complex>

#include "Parameters.h"
#include "JacobiEigenvalue.h"
#include "ToolFunctions.h"
#include "Model.h"

using namespace std;

int main(){
    
    ofstream outfile;
    outfile.open("data.txt");
    
    for(int indi=0; indi<200; indi++){
        
        //mass in MeV
        double mA=pow(10,2.+2.0/200.*indi);
        outfile<<setw(15)<<mA;
        
        double tanb=10.;
        double P11=0.03;
        double WinoMass=500000, mueff=500000;
        double lambdaSHH=0.3;    //SHuHd coupling
        
        BSModel nmssm;
        nmssm.set_model(1);
        nmssm.set_mA(mA);
        nmssm.set_tanb(tanb);
        nmssm.set_P11(P11);
        nmssm.set_WinoMass(WinoMass);
        nmssm.set_mueff(mueff);
        nmssm.set_lambdaSHH(lambdaSHH);
        
        outfile<<setw(15)<<nmssm.DecayWidth_Agaga()<<setw(15)<<nmssm.DecayWidth_ee()<<setw(15)<<nmssm.DecayWidth_mumu()<<setw(15)<<nmssm.DecayWidth_tautau()<<setw(15)<<nmssm.DecayWidth_TotalHadron()<<setw(15)<<nmssm.DecayWidth_Quark()<<setw(15)<<nmssm.DecayWidth_Gluon();
        
        outfile<<endl;
    }
    
    outfile.close();
    
    return 0;
}

