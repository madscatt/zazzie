//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_CUDAGV
#define H_CUDA_CUDAGV

#include <cuda.h>
#include <cuda_runtime.h>

#include "GV.h"
#include "cudaUtil.h"
#include "wrapperCudaKernel.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class cudaGV;
}

// definitions
class sascalc::cudaGV :
    public sascalc::GV
{
    // data
    protected:
        TYPE * _gpu_gv;
        TYPE * _gpu_q;
        TYPE * _Is_real;
        TYPE * _Is_imag;
        TYPE * _gpu_Is_real;
        TYPE * _gpu_Is_imag;

        TYPE * _gpu_coordinates;
        TYPE * _gpu_b;

    // methods
    protected:
        virtual void _updateGPU(const TYPE *const b, const int xray);
        virtual void _collectGPU();
    public:
        virtual void batch_load(const int offset, const int extend);
        virtual void calculate(const TYPE *const b, const int frame = 0, const int xray = 0);

    // meta
    public:
        cudaGV(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q);
        virtual ~cudaGV();

    // disallow copy and assign
    private:
        inline cudaGV(const cudaGV &);
        inline const cudaGV & operator=(const cudaGV &);
};

#endif
