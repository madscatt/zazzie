//
// Hailiang Zhang
// NIST & UTK
//

#include "cudaKernel_GV.h"

/// @par Main functionality
/// Calculate I(q) based on GV method
/// @param [in] gpu_coor the XYX coordinates where the leading dimension is Natoms and the slow dimension is X->Y->Z
/// @param [in] gpu_gv the gv data where the leading dimension is Ngv and the slow dimension is X->Y->Z
/// @param [in] gpu_Ireal the Ireal data where the leading dimension is Nq and the slow dimension is Ngv
__global__ void
sascalc::cuda::
cudaKernel_GV_calc(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int Nq, const int Ngv, const int xray)
{
	// get the thread id
    //int idx = blockIdx.x*blockDim.x + threadIdx.x;

    //printf("blockIdx.x: %d blockDim.x: %d threadIdx.x: %d idx: %d Natoms %d Nq %d Ngv %d\n",blockIdx.x,blockDim.x,threadIdx.x,idx,Natoms,Nq,Ngv);

    // get q values
    const int iq = threadIdx.x;
    const int igv = blockIdx.x;
    const TYPE qmag = gpu_q[iq];
    const TYPE qx = qmag * gpu_gv[igv];
    const TYPE qy = qmag * gpu_gv[Ngv + igv];
    const TYPE qz = qmag * gpu_gv[Ngv*2 + igv];

    // locals
    int iatom;
    TYPE x, y, z, b;
    TYPE q_dot_r;
    TYPE sine, cose;
    
    // summation
    TYPE Ireal = 0.0;
    TYPE Iimag = 0.0;
    for (iatom=0; iatom<Natoms; ++iatom)
    {
        x = gpu_coor[iatom];
        y = gpu_coor[Natoms + iatom];
        z = gpu_coor[Natoms*2 + iatom];
        if (xray) b = gpu_b[Natoms*iq+iatom];
        else b = gpu_b[iatom];
        q_dot_r = qx*x+qy*y+qz*z;
        sincos_TYPE(q_dot_r, &sine, &cose);
        Ireal += b*cose;
        Iimag += b*sine;
        //if (iq==1 && igv==1) printf("%8.3f %8.3f %8.3f  %8.3f  %8.3f %8.3f %8.3f\n",x,y,z,b,qx,qy,qz);
    }

    // save
    gpu_I[igv*Nq + iq] = Ireal*Ireal + Iimag*Iimag;
}

/// @par Main functionality
/// Calculate I(q) based on GV method
/// @param [in] gpu_coor the XYX coordinates where the leading dimension is Natoms and the slow dimension is X->Y->Z
/// @param [in] gpu_gv the gv data where the leading dimension is Ngv and the slow dimension is X->Y->Z
/// @param [in] gpu_Ireal the Ireal data where the leading dimension is Nq and the slow dimension is Ngv
__global__ void
sascalc::cuda::
cudaKernel_GV_calc_complex(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I_real, TYPE * const gpu_I_imag, const int Natoms, const int Nq, const int Ngv, const int xray)
{
	// get the thread id
    //int idx = blockIdx.x*blockDim.x + threadIdx.x;

    //printf("blockIdx.x: %d blockDim.x: %d threadIdx.x: %d idx: %d Natoms %d Nq %d Ngv %d\n",blockIdx.x,blockDim.x,threadIdx.x,idx,Natoms,Nq,Ngv);

    // get q values
    const int iq = threadIdx.x;
    const int igv = blockIdx.x;
    const TYPE qmag = gpu_q[iq];
    const TYPE qx = qmag * gpu_gv[igv];
    const TYPE qy = qmag * gpu_gv[Ngv + igv];
    const TYPE qz = qmag * gpu_gv[Ngv*2 + igv];

    // locals
    int iatom;
    TYPE x, y, z, b;
    TYPE q_dot_r;
    TYPE sine, cose;
    
    // summation
    TYPE Ireal = 0.0;
    TYPE Iimag = 0.0;
    for (iatom=0; iatom<Natoms; ++iatom)
    {
        x = gpu_coor[iatom];
        y = gpu_coor[Natoms + iatom];
        z = gpu_coor[Natoms*2 + iatom];
        if (xray) b = gpu_b[Natoms*iq+iatom];
        else b = gpu_b[iatom];
        q_dot_r = qx*x+qy*y+qz*z;
        sincos_TYPE(q_dot_r, &sine, &cose);
        Ireal += b*cose;
        Iimag += b*sine;
        //if (iq==1 && igv==1) printf("%8.3f %8.3f %8.3f  %8.3f  %8.3f %8.3f %8.3f\n",x,y,z,b,qx,qy,qz);
    }

    // save
    gpu_I_real[igv*Nq + iq] = Ireal*Ireal;
    gpu_I_imag[igv*Nq + iq] = Iimag*Iimag;
}

/*
/// @par Main functionality
/// Calculate I(q) based on GV method for solvent contrast
/// @param [in] gpu_coor the XYX coordinates where the leading dimension is Natoms and the slow dimension is X->Y->Z
/// @param [in] gpu_gv the gv data where the leading dimension is Ngv and the slow dimension is X->Y->Z
/// @param [in] gpu_Ireal the Ireal data where the leading dimension is Nq and the slow dimension is Ngv
__global__ void
sascalc::cuda::
cudaKernel_GV_calc_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int Nq, const int Ngv, const TYPE sld_solvent)
{
	// get the thread id
    //int idx = blockIdx.x*blockDim.x + threadIdx.x;

    //printf("blockIdx.x: %d blockDim.x: %d threadIdx.x: %d idx: %d Natoms %d Nq %d Ngv %d\n",blockIdx.x,blockDim.x,threadIdx.x,idx,Natoms,Nq,Ngv);

    // get q values
    const int iq = threadIdx.x;
    const int igv = blockIdx.x;
    const TYPE qmag = gpu_q[iq];
    const TYPE qx = qmag * gpu_gv[igv];
    const TYPE qy = qmag * gpu_gv[Ngv + igv];
    const TYPE qz = qmag * gpu_gv[Ngv*2 + igv];

    // locals
    int iatom;
    TYPE x, y, z, b, radius;
    TYPE q_dot_r;
    TYPE sine, cose;

    // summation
    TYPE Ireal = 0.0;
    TYPE Iimag = 0.0;
    for (iatom=0; iatom<Natoms; ++iatom)
    {
        x = gpu_coor[iatom];
        y = gpu_coor[Natoms + iatom];
        z = gpu_coor[Natoms*2 + iatom];
        radius = gpu_radius[iatom];
        b = 4./3.*M_PI*radius*radius*radius;
        b *= sld_solvent;
        q_dot_r = qx*x+qy*y+qz*z;
        sincos_TYPE(q_dot_r, &sine, &cose);
        Ireal += b*cose;
        Iimag += b*sine;
    }

    // save
    gpu_I[igv*Nq + iq] = Ireal*Ireal + Iimag*Iimag;
}

/// @par Main functionality
/// Calculate I(q) based on GV method for solvent contrast
/// @param [in] gpu_coor the XYX coordinates where the leading dimension is Natoms and the slow dimension is X->Y->Z
/// @param [in] gpu_gv the gv data where the leading dimension is Ngv and the slow dimension is X->Y->Z
/// @param [in] gpu_Ireal the Ireal data where the leading dimension is Nq and the slow dimension is Ngv
__global__ void
sascalc::cuda::
cudaKernel_GV_calc_Merge_GVVV(const TYPE * const gpu_coor, const TYPE * const gpu_b, const TYPE * const gpu_radius, const TYPE * const gpu_gv, const TYPE * const gpu_q, TYPE * const gpu_I, const int Natoms, const int Nq, const int Ngv, const TYPE sld_solvent, const TYPE scale_volume)
{
	// get the thread id
    //int idx = blockIdx.x*blockDim.x + threadIdx.x;

    //printf("blockIdx.x: %d blockDim.x: %d threadIdx.x: %d idx: %d Natoms %d Nq %d Ngv %d\n",blockIdx.x,blockDim.x,threadIdx.x,idx,Natoms,Nq,Ngv);

    // get q values
    const int iq = threadIdx.x;
    const int igv = blockIdx.x;
    const TYPE qmag = gpu_q[iq];
    const TYPE qx = qmag * gpu_gv[igv];
    const TYPE qy = qmag * gpu_gv[Ngv + igv];
    const TYPE qz = qmag * gpu_gv[Ngv*2 + igv];

    // locals
    int iatom;
    TYPE x, y, z, b, radius;
    TYPE q_dot_r;
    TYPE sine, cose;
    
    //const TYPE scale = 0.7; ///< @note ZHL hacked

    // summation
    TYPE Ireal = 0.0;
    TYPE Iimag = 0.0;
    for (iatom=0; iatom<Natoms; ++iatom)
    {
        x = gpu_coor[iatom];
        y = gpu_coor[Natoms + iatom];
        z = gpu_coor[Natoms*2 + iatom];
        radius = gpu_radius[iatom];
        b = gpu_b[iatom];
        //b -= 4./3.*M_PI*radius*radius*radius*sld_solvent;
        b -= scale_volume*4./3.*M_PI*radius*radius*radius*sld_solvent;
        q_dot_r = qx*x+qy*y+qz*z;
        sincos_TYPE(q_dot_r, &sine, &cose);
        Ireal += b*cose;
        Iimag += b*sine;
    }

    // save
    gpu_I[igv*Nq + iq] = Ireal*Ireal + Iimag*Iimag;
}

*/
