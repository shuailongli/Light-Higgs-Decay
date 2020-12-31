#ifndef Parameters_H
#define Parameters_H

#include <cmath>

#define PI 3.141592654

//All masses are in units of GeV
#define vev 246.0
#define fpi 0.093  //pion decay constant
#define mpi 0.135
#define meta 0.548
#define metap 0.958
#define mpi8 0.576
#define mpi9 0.942
#define mKpm 0.494
#define mK0 0.498
#define mrho 0.783

#define Delta (-0.368*0.368)
#define Bmudf (mpi*mpi) //+mKpm*mKpm-mK0*mK0) //B*mu/fpi
#define Bmddf (mpi*mpi) //-mKpm*mKpm+mK0*mK0) //B*md/fpi
#define Bmsdf (0.688*0.688)   //B*ms/fpi
#define invCf2 (0.692*0.692)  // 1/(C*fpi^2)
#define coseta 0.97437      //eta = eta8 * coseta -eta9 * sineta
#define sineta (-0.224951)  //eta' = eta8 * sineta + eta9 * coseta

#define mu 0.0022
#define md 0.0047
#define ms 0.096
#define mc 1.280
#define mb 4.180
#define mt 173.1
#define me 0.000511
#define mmu 0.10566
#define mtau 1.7768

#define mW 80.390
#define mZ 91.190
#define mh 125.



#endif
