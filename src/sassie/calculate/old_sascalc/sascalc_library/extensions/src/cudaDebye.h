//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_CUDADEBYE
#define H_CUDA_CUDADEBYE

#include <cuda.h>
#include <cuda_runtime.h>

#include "Debye.h"
#include "cudaUtil.h"
#include "wrapperCudaKernel.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class cudaDebye;
}

// definitions
class sascalc::cudaDebye :
    public sascalc::Debye
{
    // data
    protected:
        TYPE * _gpu_q;
        TYPE * _Ia;
        TYPE * _gpu_Ia;

        TYPE * _gpu_coordinates;
        TYPE * _gpu_b;

    // methods
    protected:
        virtual void _allocateGPU(const TYPE *const b, const int xray);
        virtual void _collectGPU();
    public:
        virtual void batch_load(const int offset, const int extend);
        virtual void calculate(const TYPE *const b, const int frame = 0, const int xray = 0);

    // meta
    public:
        cudaDebye(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q);
        virtual ~cudaDebye();

    // disallow copy and assign
    private:
        inline cudaDebye(const cudaDebye &);
        inline const cudaDebye & operator=(const cudaDebye &);
};

#endif
