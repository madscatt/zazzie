//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_WAPPERCUDAKERNEL
#define H_WAPPERCUDAKERNEL

#include "cudaUtil.h"
#include "cudaKernel_overlap.h"

void
wrapper_cudaOverlap(const double * const gpu_coor, const double * const gpu_r, const int natoms, const double scale, int * const gpu_noverlap, const int * const gpu_bond_table=0);


#endif
