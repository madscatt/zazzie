//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaDebye.h"
#include <math.h>
#include <iostream>
using namespace std;

//////////////////////////////////////////////////////////
/// construction
//////////////////////////////////////////////////////////
sascalc::cudaDebye::
cudaDebye(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q):
Debye(Natoms,Nframes,Nq,Nr,dr,coordinates,q)
{
    CALL_CUDA(cudaMalloc((void**)&_gpu_q, _Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_coordinates, _Natoms*3*_Nframes*sizeof(TYPE)));
    //CALL_CUDA(cudaMemcpy(_gpu_coordinates, _coordinates, _Natoms*3*_Nframes*sizeof(TYPE), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(_gpu_q, _q, _Nq*sizeof(TYPE), cudaMemcpyHostToDevice));
}

//////////////////////////////////////////////////////////
/// batch load
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
batch_load(const int offset, const int extend)
{
    _coordinates += offset*_Natoms*3;
    CALL_CUDA(cudaMemcpy(_gpu_coordinates, _coordinates, _Natoms*3*extend*sizeof(TYPE), cudaMemcpyHostToDevice));
}

//////////////////////////////////////////////////////////
/// set up for Debye
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
_allocateGPU(const TYPE *const b, const int xray)
{
    CALL_CUDA(cudaMalloc((void**)&_gpu_Ia, _Natoms*_Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMallocHost((void**)&_Ia, _Natoms*_Nq*sizeof(TYPE)));

    if (xray)
    {
        CALL_CUDA(cudaMalloc((void**)&_gpu_b, _Nq*_Natoms*sizeof(TYPE)));
        CALL_CUDA(cudaMemcpy(_gpu_b, b, _Nq*_Natoms*sizeof(TYPE), cudaMemcpyHostToDevice));
    }
    else
    {
        CALL_CUDA(cudaMalloc((void**)&_gpu_b, _Natoms*sizeof(TYPE)));
        CALL_CUDA(cudaMemcpy(_gpu_b, b, _Natoms*sizeof(TYPE), cudaMemcpyHostToDevice));
    }
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
calculate(const TYPE *const B, const int frame, const int xray)
{
    if (xray)
    {
        std::cerr<<"Xray not supported in Debye calculation yet"<<std::endl;
        exit(0);
    }
    _allocateGPU(B,xray);
    sascalc::cuda::wrapper_cudaDebye_calcIq(&_gpu_coordinates[frame*_Natoms*3], _gpu_b, _gpu_q, _gpu_Ia, _Natoms, _Nq, xray);
    _collectGPU();
}

//////////////////////////////////////////////////////////
/// harvest
////////////////////////////////////////////////////////////
void
sascalc::cudaDebye::
_collectGPU()
{
    CALL_CUDA(cudaMemcpy(_Ia, _gpu_Ia, _Natoms*_Nq*sizeof(TYPE), cudaMemcpyDeviceToHost));
    for (int i=0; i<_Nq; ++i)
    {
        _Iq[i] = 0.0;
        for (int j=0; j<_Natoms; ++j)
        {
            _Iq[i] += _Ia[i*_Natoms+j];
        }
    }
}

//////////////////////////////////////////////////////////
/// destruction
//////////////////////////////////////////////////////////
sascalc::cudaDebye::
~cudaDebye()
{
    CALL_CUDA(cudaFree(_gpu_q));
    CALL_CUDA(cudaFree(_gpu_Ia));
    CALL_CUDA(cudaFreeHost(_Ia));
    CALL_CUDA(cudaFree(_gpu_coordinates));
    CALL_CUDA(cudaFree(_gpu_b));
}
