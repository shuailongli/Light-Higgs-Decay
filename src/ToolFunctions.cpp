#include <cmath>
#include <complex>

#include "Parameters.h"
#include "ToolFunctions.h"

using namespace std;

//inpout: tau=4*mi^2/mH^2
//output: function F defined as Eq. (2.47) in hep-ph/0503172
complex<double> loop_func(double tau){
    double rea, img;
    if(tau>=1){
        rea=2*tau*asin(1/sqrt(tau))*asin(1/sqrt(tau));
        img=0;
        return complex<double>(rea,img);
    }
    else if(tau>0){
        complex<double> temp(log((1+sqrt(1-tau))/(1-sqrt(1-tau))),-PI);
        return complex<double>(-0.5*tau,0)*temp*temp;
    }
    else return 0;
}

//input: three 3 x 3 matrix mat1, mat2, mat3
//output: mat1 * mat2 * mat3 + (all permutations of {1, 2, 3}), an 3 x 3 matrix
// "*" denotes the matrix product operator
void PermuProduct(double mat1[3][3], double mat2[3][3], double mat3[3][3], double remat[3][3]){
    
    double Matrix[3][3][3];
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            Matrix[0][i][j]=mat1[i][j];
            Matrix[1][i][j]=mat2[i][j];
            Matrix[2][i][j]=mat3[i][j];
        }
    }
    
    //define a matrix storing 6 permutations of {1,2,3}
    int index[6][3]={{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++) remat[i][j]=0;
    }
    
    for(int i=0;i<6;i++){
        for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
                for(int i1=0;i1<3;i1++){
                    for(int i2=0;i2<3;i2++){
                        remat[j][k]+=Matrix[index[i][0]][j][i1] *Matrix[index[i][1]][i1][i2] *Matrix[index[i][2]][i2][k];
                    }
                }
            }
        }
    }
}


//input: four 3 x 3 matrix mat1, mat2, mat3, mat4
//output: mat1 * mat2 * mat3 * mat4 + (all permutations of {1, 2, 3, 4}), an 3 x 3 matrix
// "*" denotes the matrix product operator
void PermuProduct(double mat1[3][3], double mat2[3][3], double mat3[3][3], double mat4[3][3], double remat[3][3]){
    
    double Matrix[4][3][3];
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            Matrix[0][i][j]=mat1[i][j];
            Matrix[1][i][j]=mat2[i][j];
            Matrix[2][i][j]=mat3[i][j];
            Matrix[3][i][j]=mat4[i][j];
        }
    }
    
    //define a matrix storing 6 permutations of {1,2,3}
    int index[24][4]={{0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1,
        2}, {0, 3, 2, 1}, {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2,
            3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0}, {2, 0, 1, 3}, {2, 0, 3, 1}, {2,
                1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0}, {3, 0, 1,
                    2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2,
                        1, 0}};
    
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++) remat[i][j]=0;
    }
    
    for(int i=0;i<24;i++){
        for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
                for(int i1=0;i1<3;i1++){
                    for(int i2=0;i2<3;i2++){
                        for(int i3=0;i3<3;i3++){
                            remat[j][k]+=Matrix[index[i][0]][j][i1] *Matrix[index[i][1]][i1][i2] *Matrix[index[i][2]][i2][i3] *Matrix[index[i][3]][i3][k];
                        }
                    }
                }
            }
        }
    }
}



//numerical integration of eq. (22) in 1612.06538
double appp_func(double a, double mP, double m1, double m2, double m3){
    return sqrt(abs(1-2*(m2*m2+m3*m3)/a+(m2*m2-m3*m3)*(m2*m2-m3*m3)/(a*a))) *sqrt(abs( (1+(a-m1*m1)/(mP*mP))*(1+(a-m1*m1)/(mP*mP))-(4*a)/(mP*mP)) );
}
double APPP(double mP, double m1, double m2, double m3){
    if(mP<=(m1+m2+m3)) return 0;
    else{
        double a=(m2+m3)*(m2+m3), b=(mP-m1)*(mP-m1);
        double dx=(b-a)/10000;
        double sum=0;
        for(double x = a; x < b; x += dx)
        {
            sum += appp_func(x,mP,m1,m2,m3) + 4*appp_func(x + dx/2,mP,m1,m2,m3) + appp_func(x + dx,mP,m1,m2,m3);
        }
        
        sum *= dx/6;
        
        return sum;
    }
}



//input: mA, invariant mass square s;
//output: Beta, Betap: anomaly scalar
pair <complex<double>, complex<double> > BEta(double mA, double s){
    double Fpid8=1/1.5323389627849169, Fpid0 = 1/1.1584711549017292;
    double thetaP=thetaP = -20*PI/180.;
    
    complex<double> Form((1 - s/(4*mpi*mpi))* sqrt((s - 4*mpi*mpi)/s)* log((1 + sqrt((s - 4*mpi*mpi)/s))/(1 - sqrt((s - 4*mpi*mpi)/s))) - 2, (1 - s/(4*mpi*mpi))* sqrt((s - 4*mpi*mpi)/s)*PI);
    complex<double> D1=1 - s/(mrho*mrho) - s/(96 *PI*PI *fpi*fpi)* log(mrho*mrho/(mpi*mpi)) - mpi*mpi/(24 *PI*PI *fpi*fpi)* Form;
    
    complex<double> Bppgeta=3/(12*sqrt(3.) *PI*PI *fpi*fpi*fpi) *(Fpid8 *cos(thetaP) -sqrt(2.)* Fpid0 *sin(thetaP)) *((1 + s/(2* mrho*mrho))/D1);
    complex<double> Bppgetap=3/(12*sqrt(3.) *PI*PI *fpi*fpi*fpi) *(Fpid8 *sin(thetaP) +sqrt(2.)* Fpid0 *cos(thetaP)) *((1 + s/(2* mrho*mrho))/D1);
    
    return make_pair(Bppgeta, Bppgetap);
}

double Gamma0(double mA, double s){
    return (0.302862*0.302862 *s *sqrt((s - 4 *mpi*mpi)/s) *sqrt((s - 4 *mpi*mpi)/s) *sqrt((s - 4 *mpi*mpi)/s) *(mA*mA - s) *(mA*mA - s) *(mA*mA - s))/(12 *(8 *PI *mA) *(8 *PI *mA) *(8 *PI *mA));
}


//return the sign of x
int SIGN(double x){
    if(x>0) return 1;
    else if(x==0) return 0;
    else return -1;
}
