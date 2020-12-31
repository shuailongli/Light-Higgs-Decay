#ifndef ToolFunctions_H
#define ToolFunctions_H

#include <complex>

using namespace std;

//inpout: tau=4*mi^2/mH^2
//output: F function defined as Eq. (2.47) in hep-ph/0503172
complex<double> loop_func(double tau);

//input: three 3 x 3 matrix mat1, mat2, mat3
//output: mat1 * mat2 * mat3 + (all permutations of {1, 2, 3}), an 3 x 3 matrix
// "*" denotes the matrix product operator
void PermuProduct(double mat1[3][3], double mat2[3][3], double mat3[3][3], double remat[3][3]);

//input: four 3 x 3 matrix mat1, mat2, mat3, mat4
//output: mat1 * mat2 * mat3 * mat4 + (all permutations of {1, 2, 3, 4}), an 3 x 3 matrix
// "*" denotes the matrix product operator
void PermuProduct(double mat1[3][3], double mat2[3][3], double mat3[3][3], double mat4[3][3], double remat[3][3]);

//numerical integration of eq. (22) in 1612.06538
double APPP(double mP, double m1, double m2, double m3);


//input: mA, invariant mass square s;
//output: Beta, Betap: anomaly scalar
pair <complex<double>, complex<double> > BEta(double mA, double s);
double Gamma0(double mA, double s);


//return the sign of x
int SIGN(double x);


double abs_Delta_K(double x);
double abs_Delta_Pi(double x);
double abs_Gamma_K(double x);
double abs_Gamma_Pi(double x);
double abs_Theta_K(double x);
double abs_Theta_Pi(double x);
double alphasQ(double x);
double arg_Delta_K(double x);
double arg_Delta_Pi(double x);
double arg_Gamma_K(double x);
double arg_Gamma_Pi(double x);
double arg_Theta_K(double x);
double arg_Theta_Pi(double x);

#endif
