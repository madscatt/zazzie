//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaDebyePBC.h"
#include <math.h>
#include <iostream>
using namespace std;

//////////////////////////////////////////////////////////
/// construction
//////////////////////////////////////////////////////////
sascalc::cudaDebyePBC::
cudaDebyePBC(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q, const TYPE box):
cudaDebye(Natoms,Nframes,Nq,Nr,dr,coordinates,q),
_boxl(box)
{
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaDebyePBC::
calculate(const TYPE *const B, const int frame, const int xray)
{
    if (xray)
    {
        std::cerr<<"Xray not supported in DebyePBC calculation yet"<<std::endl;
        exit(0);
    }
    _allocateGPU(B,xray);
    sascalc::cuda::wrapper_cudaDebyePBC_calcIq(_boxl, &_gpu_coordinates[frame*_Natoms*3], _gpu_b, _gpu_q, _gpu_Ia, _Natoms, _Nq, xray);
    _collectGPU();
}

//////////////////////////////////////////////////////////
/// destruction
//////////////////////////////////////////////////////////
sascalc::cudaDebyePBC::
~cudaDebyePBC()
{
}
