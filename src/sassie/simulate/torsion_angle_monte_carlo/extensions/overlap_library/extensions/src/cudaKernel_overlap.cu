//
// Hailiang Zhang
// NIST & UTK
//

#include <cudaKernel_overlap.h>

/// @par Main functionality
/// check for overlap
/// @param [in] gpu_coor the XYX coordinates where the leading dimension is X->Y->Z and slow dimension is natoms
__global__ void cudaKernel_overlap(const double * const gpu_coor, const double * const gpu_radius, const int natoms, const int npairs, const double scale, int * const gpu_noverlap, const int * const gpu_bond_table)
{
	// get the thread id
    int idx = blockIdx.x*blockDim.x + threadIdx.x;

    // thread boundary check
    if (idx>=npairs) return;

    // return if bonded
    if (gpu_bond_table && gpu_bond_table[idx]) return;

    // get the pair indices
    /*
    int idx_pair_small, idx_pair_big=0;
    int npairs_sofar = 0;
    do
    {
        ++idx_pair_big;
        npairs_sofar += idx_pair_big;
    }while (npairs_sofar < idx+1);
    npairs_sofar -= idx_pair_big;
    idx_pair_small = idx-npairs_sofar;
    */
    int idx_pair_big = ceil((1.+sqrt(1.+8.*(idx+1)))/2.)-1;
    int idx_pair_small = idx - ((idx_pair_big-1)*idx_pair_big)/2;
    //if (idx<50) printf("idx: %d idx_pair_big: %d idx_pair_small: %d\n",idx,idx_pair_big,idx_pair_small);

    // get the pair coordinates and radii
    double xij = gpu_coor[3*idx_pair_big] - gpu_coor[3*idx_pair_small];
    double yij = gpu_coor[3*idx_pair_big+1] - gpu_coor[3*idx_pair_small+1];
    double zij = gpu_coor[3*idx_pair_big+2] - gpu_coor[3*idx_pair_small+2];
    double ri = gpu_radius[idx_pair_big];
    double rj = gpu_radius[idx_pair_small];

    // get the distance
    double r = sqrt(xij*xij + yij*yij + zij*zij);
    //if (r<(ri+rj)*scale) printf("Overlap found. idx: %d, idx_pair_big: %d, idx_pair_small: %d, ri: %f, rj: %f, r: %f, scale: %f\n",idx,idx_pair_big,idx_pair_small,ri,rj,r,scale);
    if (r<(ri+rj)*scale) atomicAdd(gpu_noverlap,1);

    // return
    return;
}
