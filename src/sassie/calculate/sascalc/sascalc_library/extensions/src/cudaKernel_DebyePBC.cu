//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaKernel_DebyePBC.h"

/// @par Main functionality
/// Calculate I(q) based on DebyePBC summation
/// @param [in] gpu_coor the XYX coordinates where the leading dimension is Natoms and the slow dimension is X->Y->Z
/// @param [in] gpu_Ia the partial intensity where the leading dimension is Natoms and the slow dimension is Nq
__global__ void
sascalc::cuda::
cudaKernel_DebyePBC_calc(const TYPE boxl, const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_q, TYPE * const gpu_Ia, const int Natoms, const int Nq, const int Xray)
{
	// get the thread id and return if necessary
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx>=Natoms) return;

    // shared memory
    __shared__ TYPE shared_x[BLOCKDIM];
    __shared__ TYPE shared_y[BLOCKDIM];
    __shared__ TYPE shared_z[BLOCKDIM];
    __shared__ TYPE shared_b[BLOCKDIM];

    // coordinate/b for this thread
    const TYPE xi = gpu_coor[idx];
    const TYPE yi = gpu_coor[Natoms+idx];
    const TYPE zi = gpu_coor[Natoms*2+idx];
    const TYPE bi = gpu_b[idx];

    // local variables
    int idx_atom,idx_block, idx_q;
    TYPE r,q,qr,bj;

    for (idx_q=0; idx_q<Nq; ++idx_q) gpu_Ia[idx_q*Natoms + idx] = bi*bi;

    // loop over number of blocks
    //for (idx_block=blockIdx.x; idx_block<((Natoms-1)/blockDim.x)+1; ++idx_block)
    for (idx_block=blockIdx.x; idx_block<gridDim.x; ++idx_block)
    {
        // load coordinates and bs into shared memory
        __syncthreads();
        idx_atom = idx_block*blockDim.x + threadIdx.x;
        shared_x[threadIdx.x] = gpu_coor[idx_atom];
        shared_y[threadIdx.x] = gpu_coor[Natoms + idx_atom];
        shared_z[threadIdx.x] = gpu_coor[Natoms*2 + idx_atom];
        shared_b[threadIdx.x] = gpu_b[idx_atom];
        __syncthreads();
        for (idx_atom=0; idx_atom<blockDim.x; ++idx_atom)
        {
            if (idx_block==blockIdx.x && idx_atom<=threadIdx.x) continue;
            if (idx_block==gridDim.x-1) {if (idx_block*blockDim.x + idx_atom >= Natoms) break;}
            r = sqrt ( powf(xi-shared_x[idx_atom], 2.0) + powf(yi-shared_y[idx_atom], 2.0) + powf(zi-shared_z[idx_atom], 2.0) );
            r -= boxl*(round(r/boxl));
            bj = shared_b[idx_atom];
            for (idx_q=0; idx_q<Nq; ++idx_q)
            {
                q = gpu_q[idx_q];
                qr = q*r;
                if (qr==0.0) gpu_Ia[idx_q*Natoms + idx] += 2.*bi*bj;
                else gpu_Ia[idx_q*Natoms + idx] += 2.*bi*bj*sinf(q*r)/(q*r);
                //if (idx>=1958 && idx_q==0) printf("self: %04d, other: %04d, qr: %f, bi: %f, bj: %f, gpu_Ia[idx_q*Natoms + idx]: %f\n",idx,idx_block*blockDim.x + idx_atom,qr,bi,bj,gpu_Ia[idx_q*Natoms + idx]);
            }
         }
    }
}

