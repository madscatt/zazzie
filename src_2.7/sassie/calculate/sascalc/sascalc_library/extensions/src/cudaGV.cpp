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
cudaGV(const int Natoms, const ScVars & scvar, const double *const B_neutron_array, const double *const B_xray_array):
GV(Natoms, scvar, B_neutron_array, B_xray_array)
{
    const int Nitems_neutron = _B_neutron_vector.size();
    const int Nitems_xray = _B_xray_vector.size();
    CALL_CUDA(cudaMalloc((void**)&_gpu_q, _Nq*sizeof(double)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_coor, 3*_Natoms*sizeof(double)));
    CALL_CUDA(cudaMemcpy(_gpu_q, _q, _Nq*sizeof(double), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMalloc((void**)&_gpu_B_neutron, Nitems_neutron*_Natoms*sizeof(double)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_B_xray, Nitems_xray*_Natoms*_Nq*sizeof(double)));
    CALL_CUDA(cudaMemcpy(_gpu_B_neutron, B_neutron_array, Nitems_neutron*_Natoms*sizeof(double), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(_gpu_B_xray, B_xray_array, Nitems_xray*_Natoms*_Nq*sizeof(double), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMalloc((void**)&_gpu_flag_skip, (Nitems_neutron+Nitems_xray)*sizeof(int)));
}

//////////////////////////////////////////////////////////
/// set up 
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
_updateGPU(const int *const flag_skip, const int Ngv, const double *const coor)
{
    const int Nitems = _B_neutron_vector.size()+_B_xray_vector.size();
    CALL_CUDA(cudaMallocHost((void**)&_Is, Nitems*Ngv*_Nq*sizeof(double)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_Is, Nitems*Ngv*_Nq*sizeof(double)));
    CALL_CUDA(cudaMalloc((void**)&_gpu_gv, 3*Ngv*sizeof(double)));
    double *gv = _getGV(Ngv);
    CALL_CUDA(cudaMemcpy(_gpu_gv, gv, 3*Ngv*sizeof(double), cudaMemcpyHostToDevice));
    delete [] gv;
    CALL_CUDA(cudaMemcpy(_gpu_coor, coor, 3*_Natoms*sizeof(double), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(_gpu_flag_skip, flag_skip, Nitems*sizeof(int), cudaMemcpyHostToDevice));
}

//////////////////////////////////////////////////////////
/// calcIq 
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
_calculate_singleFrame_fixed(const int *const flag_skip, const int Nitems, const double *const coor, double *const Iq, const int Ngv)
{
    const int Ncontrast_neutron = _B_neutron_vector.size()/3;
    const int Ncontrast_xray = _B_xray_vector.size()/3;
    _updateGPU(flag_skip, Ngv, coor);
    sascalc::cuda::wrapper_cudaGV_calcIq(_gpu_coor, _gpu_B_neutron, _gpu_B_xray, _gpu_flag_skip, _gpu_gv, _gpu_q, _gpu_Is, _Natoms, _Nq, Ngv, Ncontrast_neutron, Ncontrast_xray);
    _collectGPU(Ngv, flag_skip, Iq);
}

//////////////////////////////////////////////////////////
/// harvest
//////////////////////////////////////////////////////////
void
sascalc::cudaGV::
_collectGPU(const int Ngv, const int *const flag_skip, double *const Iq)
{
    const int Nitems = _B_neutron_vector.size()+_B_xray_vector.size();
    const int Ncontrast_neutron = _B_neutron_vector.size()/3;
    const int Ncontrast_xray = _B_xray_vector.size()/3;
    CALL_CUDA(cudaMemcpy(_Is, _gpu_Is, Nitems*Ngv*_Nq*sizeof(double), cudaMemcpyDeviceToHost));
    // summation
    int i,j;
    double *pI;
    int item = 0, offset=0;
    int idx_contrast, idx;
    for (idx_contrast=0; idx_contrast!=Ncontrast_neutron; ++idx_contrast)
    {
        for (idx=0; idx!=3; ++idx)
        {
            pI = &Iq[item*_Nq]; 
            if (!flag_skip[item])
            {
                for (i=0; i<_Nq; ++i)
                {
                    pI[i] = 0.0;
                    for (j=0; j<Ngv; ++j)
                    {
                        pI[i] += _Is[offset+j*_Nq+i];
                    }
                    pI[i] /= Ngv;
                }
            }
            offset += _Nq*Ngv;
            ++item;
        }
    }
    for (idx_contrast=0; idx_contrast!=Ncontrast_xray; ++idx_contrast)
    {
        for (idx=0; idx!=3; ++idx)
        {
            pI = &Iq[item*_Nq]; 
            if (!flag_skip[item])
            {
                for (i=0; i<_Nq; ++i)
                {
                    pI[i] = 0.0;
                    for (j=0; j<Ngv; ++j)
                    {
                        pI[i] += _Is[offset+j*_Nq+i];
                    }
                    pI[i] /= Ngv;
                }
            }
            offset += _Nq*Ngv;
            ++item;
        }
    }

    // free temp
    CALL_CUDA(cudaFree(_gpu_Is));
    CALL_CUDA(cudaFree(_gpu_gv));
    CALL_CUDA(cudaFreeHost(_Is));
}

//////////////////////////////////////////////////////////
/// destruction
//////////////////////////////////////////////////////////
sascalc::cudaGV::
~cudaGV()
{
    CALL_CUDA(cudaFree(_gpu_q));
    CALL_CUDA(cudaFree(_gpu_coor));
    CALL_CUDA(cudaFree(_gpu_B_neutron));
    CALL_CUDA(cudaFree(_gpu_B_xray));
    CALL_CUDA(cudaFree(_gpu_flag_skip));
    CALL_CUDA(cudaDeviceReset());
}
