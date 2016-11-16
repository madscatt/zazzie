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
cudaKernel_GV_calc(const double * const gpu_coor, const double * const gpu_B_neutron, const double * const gpu_B_xray, const int *const gpu_flag_skip, const double * const gpu_gv, const double * const gpu_q, double * const gpu_Is, const int Natoms, const int Nq, const int Ngv, const int Ncontrast_neutron, const int Ncontrast_xray)
{
	// get the thread id
    //int idx = blockIdx.x*blockDim.x + threadIdx.x;

    //printf("blockIdx.x: %d blockDim.x: %d threadIdx.x: %d idx: %d Natoms %d Nq %d Ngv %d\n",blockIdx.x,blockDim.x,threadIdx.x,idx,Natoms,Nq,Ngv);

    // get q values
    const int iq = threadIdx.x;
    const int igv = blockIdx.x;
    const double qmag = gpu_q[iq];
    const double qx = qmag * gpu_gv[igv];
    const double qy = qmag * gpu_gv[Ngv + igv];
    const double qz = qmag * gpu_gv[Ngv*2 + igv];

    // locals
    int iatom;
    double x, y, z, b;
    double q_dot_r;
    double sine, cose;
    int idx, idx_contrast;
    double Ireal;
    double Iimag;
    
    // summation
    int item = 0;
    int offset = 0;
    for (idx_contrast=0; idx_contrast!=Ncontrast_neutron; ++idx_contrast)
    {
        for (idx=0; idx!=3; ++idx)
        {
            if (!gpu_flag_skip[item])
            {
                Ireal = 0.0;
                Iimag = 0.0;
                for (iatom=0; iatom<Natoms; ++iatom)
                {
                    x = gpu_coor[iatom];
                    y = gpu_coor[Natoms + iatom];
                    z = gpu_coor[Natoms*2 + iatom];
                    b = gpu_B_neutron[idx*(Ncontrast_neutron*Natoms) + idx_contrast*Natoms + iatom];
                    q_dot_r = qx*x+qy*y+qz*z;
                    sincos(q_dot_r, &sine, &cose);
                    Ireal += b*cose;
                    Iimag += b*sine;
                    //if (iq==1 && igv==1 && idx_contrast==0 && idx==0 && iatom==1) printf("%8.3f %8.3f %8.3f  %8.3f  %8.3f %8.3f %8.3f\n",x,y,z,b,qx,qy,qz);
                }
                // save
                gpu_Is[offset + igv*Nq + iq] = Ireal*Ireal + Iimag*Iimag;
            }
            offset += Nq*Ngv;
            ++item;
        }
    }
    for (idx_contrast=0; idx_contrast!=Ncontrast_xray; ++idx_contrast)
    {
        for (idx=0; idx!=3; ++idx)
        {
            if (!gpu_flag_skip[item])
            {
                Ireal = 0.0;
                Iimag = 0.0;
                for (iatom=0; iatom<Natoms; ++iatom)
                {
                    x = gpu_coor[iatom];
                    y = gpu_coor[Natoms + iatom];
                    z = gpu_coor[Natoms*2 + iatom];
                    b = gpu_B_xray[idx*(Ncontrast_xray*Nq*Natoms) + idx_contrast*Nq*Natoms + iq*Natoms + iatom];
                    q_dot_r = qx*x+qy*y+qz*z;
                    sincos(q_dot_r, &sine, &cose);
                    Ireal += b*cose;
                    Iimag += b*sine;
                    //if (iq==1 && igv==1) printf("%8.3f %8.3f %8.3f  %8.3f  %8.3f %8.3f %8.3f\n",x,y,z,b,qx,qy,qz);
                }
                // save
                gpu_Is[offset + igv*Nq + iq] = Ireal*Ireal + Iimag*Iimag;
            }
            offset += Nq*Ngv;
            ++item;
        }
    }
}
