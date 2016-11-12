//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaGV.h"
#include <math.h>
#include <iostream>
using namespace std;

//////////////////////////////////////////////////////////
/// construction
//////////////////////////////////////////////////////////
sascalc::cudaGV::
cudaGV(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q):
GV(Natoms,Nframes,Nq,Nr,dr,coordinates,q)
{
    _Ngv=1;
    CALL_CUDA(cudaMalloc((void**)&_gpu_q, _Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_coordinates, _Natoms*3*_Nframes*sizeof(TYPE)));
    //CALL_CUDA(cudaMemcpy(_gpu_coordinates, _coordinates, _Natoms*3*_Nframes*sizeof(TYPE), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(_gpu_q, _q, _Nq*sizeof(TYPE), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Is_real, _Ngv*_Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMallocHost((void**)&_Is_real, _Ngv*_Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Is_imag, _Ngv*_Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMallocHost((void**)&_Is_imag, _Ngv*_Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_b, _Nq*sizeof(TYPE)));

    CALL_CUDA(cudaMalloc((void**)&_gpu_gv, 3*_Ngv*sizeof(TYPE)));
}

//////////////////////////////////////////////////////////
/// batch load
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
batch_load(const int offset, const int extend)
{
    _coordinates += offset*_Natoms*3;
    CALL_CUDA(cudaMemcpy(_gpu_coordinates, _coordinates, _Natoms*3*extend*sizeof(TYPE), cudaMemcpyHostToDevice));
}

//////////////////////////////////////////////////////////
/// set up 
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
_updateGPU(const TYPE *const b, const int xray)
{
    CALL_CUDA(cudaFreeHost(_Is_real));
    CALL_CUDA(cudaFreeHost(_Is_imag));
    CALL_CUDA(cudaFree(_gpu_Is_real));
    CALL_CUDA(cudaFree(_gpu_Is_imag));
    CALL_CUDA(cudaFree(_gpu_gv));
    CALL_CUDA(cudaFree(_gpu_b));

    CALL_CUDA(cudaMallocHost((void**)&_Is_real, _Ngv*_Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMallocHost((void**)&_Is_imag, _Ngv*_Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Is_real, _Ngv*_Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Is_imag, _Ngv*_Nq*sizeof(TYPE)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_gv, 3*_Ngv*sizeof(TYPE)));
    CALL_CUDA(cudaMemcpy(_gpu_gv, _gv, 3*_Ngv*sizeof(TYPE), cudaMemcpyHostToDevice));

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
sascalc::cudaGV::
calculate(const TYPE *const B, const int frame, const int xray)
{
    _updateGPU(B,xray);
    sascalc::cuda::wrapper_cudaGV_calcIq_complex(&_gpu_coordinates[frame*_Natoms*3], _gpu_b, _gpu_gv, _gpu_q, _gpu_Is_real, _gpu_Is_imag, _Natoms, _Nq, _Ngv, xray);
    _collectGPU();
}

//////////////////////////////////////////////////////////
/// harvest
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
_collectGPU()
{
    CALL_CUDA(cudaMemcpy(_Is_real, _gpu_Is_real, _Ngv*_Nq*sizeof(TYPE), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(_Is_imag, _gpu_Is_imag, _Ngv*_Nq*sizeof(TYPE), cudaMemcpyDeviceToHost));
    for (int i=0; i<_Nq; ++i)
    {
        _Iq_real[i] = 0.0;
        _Iq_imag[i] = 0.0;
        for (int j=0; j<_Ngv; ++j)
        {
            _Iq_real[i] += _Is_real[j*_Nq+i];
            _Iq_imag[i] += _Is_imag[j*_Nq+i];
        }
        _Iq[i] = (_Iq_real[i]+_Iq_imag[i])/_Ngv;
        _Iq_real[i] = pow(_Iq_real[i]/_Ngv, 0.5);
        _Iq_imag[i] = pow(_Iq_imag[i]/_Ngv, 0.5);
    }
}

//////////////////////////////////////////////////////////
/// destruction
//////////////////////////////////////////////////////////
sascalc::cudaGV::
~cudaGV()
{
    CALL_CUDA(cudaFree(_gpu_gv));
    CALL_CUDA(cudaFree(_gpu_q));
    CALL_CUDA(cudaFree(_gpu_Is_real));
    CALL_CUDA(cudaFree(_gpu_Is_imag));
    CALL_CUDA(cudaFreeHost(_Is_real));
    CALL_CUDA(cudaFreeHost(_Is_imag));
    CALL_CUDA(cudaFree(_gpu_coordinates));
    CALL_CUDA(cudaFree(_gpu_b));
    CALL_CUDA(cudaDeviceReset());
}
