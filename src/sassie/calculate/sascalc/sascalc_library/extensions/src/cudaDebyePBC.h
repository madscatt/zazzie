//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_CUDADEBYEPBC
#define H_CUDA_CUDADEBYEPBC

#include <cuda.h>
#include <cuda_runtime.h>

#include "cudaDebye.h"
#include "cudaUtil.h"
#include "wrapperCudaKernel.h"

// declaration
namespace sascalc
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class cudaDebyePBC;
}

// definitions
class sascalc::cudaDebyePBC :
    public sascalc::cudaDebye
{
    // data
    protected:
        const TYPE _boxl;

    // methods
    public:
        virtual void calculate(const TYPE *const b, const int frame = 0, const int xray = 0);

    // meta
    public:
        cudaDebyePBC(const int Natoms, const int Nframes, const int Nq, const int Nr, const TYPE dr, const TYPE *const coordinates, const TYPE *const q, const TYPE boxl);
        virtual ~cudaDebyePBC();

    // disallow copy and assign
    private:
        inline cudaDebyePBC(const cudaDebyePBC &);
        inline const cudaDebyePBC & operator=(const cudaDebyePBC &);
};

#endif
