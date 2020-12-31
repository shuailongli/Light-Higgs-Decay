#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <math.h>
#include <cmath>

#include "Model.h"


using namespace std;

int main(){
    
    ofstream outfile;
    outfile.open("data.txt");
    
    for(int i=0;i<50;i++){
        double mH=pow(10,-2+3./50.*i);
        for(int j=0;j<50;j++){
            double tanb=pow(10,1.69897/50.*j);
            double cosba=1/tanb;
            double mA=1.0, mHC=600;
            double lambdav=0;
            
            BSModel model;
            model.set_model(2);
            model.set_cosba(cosba);
            model.set_tanb(tanb);
            model.set_mH(mH);
            model.set_mA(mA);
            model.set_mHC(mHC);
            model.set_lambdav(lambdav);
            
            outfile<<setw(10)<<mH<<setw(10)<<tanb;
            
            //CP-odd Higgs decay
            outfile<<setw(15)<<model.DecayWidth_Aee()<<setw(15)<<model.DecayWidth_AGluon();
            outfile<<setw(15)<<model.BranchRatio_Aee()<<setw(15)<<model.BranchRatio_AGluon();
            
            //CP-even Higgs decay
            outfile<<setw(15)<<model.DecayWidth_Hee()<<setw(15)<<model.DecayWidth_HGluon();
            outfile<<setw(15)<<model.BranchRatio_Hee()<<setw(15)<<model.BranchRatio_HGluon();
            
            outfile<<endl;
        }
    }
}
