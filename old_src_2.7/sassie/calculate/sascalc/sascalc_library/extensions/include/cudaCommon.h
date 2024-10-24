//
// Hailiang Zhang
// NIST & UTK
//

#ifndef H_CUDA_COMMON
#define H_CUDA_COMMON

#define USE_DOUBLE

#ifdef USE_DOUBLE
#define TYPE double
#define MPI_TYPE MPI_DOUBLE
#define gsl_vector_TYPE gsl_vector
#define gsl_vector_TYPE_fprintf gsl_vector_fprintf
#define gsl_vector_TYPE_alloc gsl_vector_alloc
#define gsl_vector_TYPE_calloc gsl_vector_calloc
#define gsl_vector_TYPE_memcpy gsl_vector_memcpy
#define gsl_vector_TYPE_set gsl_vector_set
#define gsl_vector_TYPE_get gsl_vector_get
#define gsl_vector_TYPE_set_zero gsl_vector_set_zero
#define gsl_vector_TYPE_set_all gsl_vector_set_all
#define gsl_vector_TYPE_sub gsl_vector_sub
#define gsl_vector_TYPE_scale gsl_vector_scale
#define gsl_vector_TYPE_free gsl_vector_free
#define gsl_vector_TYPE_view gsl_vector_view
#define gsl_vector_TYPE_fread gsl_vector_fread
#define gsl_vector_TYPE_reverse gsl_vector_reverse
#define gsl_vector_TYPE_view gsl_vector_view
#define gsl_vector_TYPE_subvector gsl_vector_subvector
#define gsl_matrix_TYPE gsl_matrix
#define gsl_matrix_TYPE_alloc gsl_matrix_alloc
#define gsl_matrix_TYPE_calloc gsl_matrix_calloc
#define gsl_matrix_TYPE_memcpy gsl_matrix_memcpy
#define gsl_matrix_TYPE_transpose_memcpy gsl_matrix_transpose_memcpy
#define gsl_matrix_TYPE_set gsl_matrix_set
#define gsl_matrix_TYPE_get gsl_matrix_get
#define gsl_matrix_TYPE_set_zero gsl_matrix_set_zero
#define gsl_matrix_TYPE_add gsl_matrix_add
#define gsl_matrix_TYPE_scale gsl_matrix_scale
#define gsl_matrix_TYPE_free gsl_matrix_free
#define gsl_matrix_TYPE_view gsl_matrix_view
#define gsl_matrix_TYPE_row gsl_matrix_row
#define gsl_matrix_TYPE_column gsl_matrix_column
#define gsl_matrix_TYPE_get_row gsl_matrix_get_row
#define gsl_matrix_TYPE_get_col gsl_matrix_get_col
#define gsl_matrix_TYPE_set_row gsl_matrix_set_row
#define gsl_matrix_TYPE_set_col gsl_matrix_set_col
#define gsl_matrix_TYPE_const_view gsl_matrix_const_view
#define gsl_matrix_TYPE_const_submatrix gsl_matrix_const_submatrix
#define gsl_matrix_TYPE_fread gsl_matrix_fread
#define gsl_matrix_TYPE_fwrite gsl_matrix_fwrite
#define gsl_matrix_TYPE_submatrix gsl_matrix_submatrix
#define gsl_matrix_TYPE_mul_elements gsl_matrix_mul_elements
#define gsl_matrix_TYPE_div_elements gsl_matrix_div_elements
#define gsl_matrix_TYPE_add_constant gsl_matrix_add_constant
#define gsl_stats_TYPE_min gsl_stats_min
#define gsl_stats_TYPE_max gsl_stats_max
#define gsl_stats_TYPE_mean gsl_stats_mean
#define gsl_stats_TYPE_wmean gsl_stats_wmean
#define gsl_stats_TYPE_sd gsl_stats_sd
#define gsl_stats_TYPE_median_from_sorted_data gsl_stats_median_from_sorted_data
#define gsl_sort_vector_TYPE gsl_sort_vector
#define gsl_sort_vector_TYPE_index gsl_sort_vector_index
#define gsl_permute_vector_TYPE gsl_permute_vector
#define gsl_blas_TYPEcopy gsl_blas_dcopy
#define gsl_blas_TYPEaxpy gsl_blas_daxpy
#define gsl_blas_TYPEsyr gsl_blas_dsyr
#define gsl_blas_TYPEtrmv gsl_blas_dtrmv
#define gsl_blas_TYPEtrmm gsl_blas_dtrmm
#define gsl_blas_TYPEgemv gsl_blas_dgemv
#define gsl_blas_TYPEgemm gsl_blas_dgemm
#define gsl_blas_TYPEsymm gsl_blas_dsymm
#define gsl_blas_TYPEdot gsl_blas_ddot
#define gsl_blas_TYPEasum gsl_blas_dasum
#define gsl_blas_TYPEnrm2 gsl_blas_dnrm2
#define NATIVE_TYPE NATIVE_DOUBLE
#else
#define TYPE float
#define MPI_TYPE MPI_FLOAT
#define sincos_TYPE sincosf
#define gsl_vector_TYPE gsl_vector_float
#define gsl_vector_TYPE_fprintf gsl_vector_float_fprintf
#define gsl_vector_TYPE_alloc gsl_vector_float_alloc
#define gsl_vector_TYPE_calloc gsl_vector_float_calloc
#define gsl_vector_TYPE_memcpy gsl_vector_float_memcpy
#define gsl_vector_TYPE_set gsl_vector_float_set
#define gsl_vector_TYPE_get gsl_vector_float_get
#define gsl_vector_TYPE_set_zero gsl_vector_float_set_zero
#define gsl_vector_TYPE_set_all gsl_vector_float_set_all
#define gsl_vector_TYPE_sub gsl_vector_float_sub
#define gsl_vector_TYPE_scale gsl_vector_float_scale
#define gsl_vector_TYPE_free gsl_vector_float_free
#define gsl_vector_TYPE_view gsl_vector_float_view
#define gsl_vector_TYPE_fread gsl_vector_float_fread
#define gsl_vector_TYPE_reverse gsl_vector_float_reverse
#define gsl_vector_TYPE_view gsl_vector_float_view
#define gsl_vector_TYPE_subvector gsl_vector_float_subvector
#define gsl_matrix_TYPE gsl_matrix_float
#define gsl_matrix_TYPE_alloc gsl_matrix_float_alloc
#define gsl_matrix_TYPE_calloc gsl_matrix_float_calloc
#define gsl_matrix_TYPE_memcpy gsl_matrix_float_memcpy
#define gsl_matrix_TYPE_transpose_memcpy gsl_matrix_float_transpose_memcpy
#define gsl_matrix_TYPE_set gsl_matrix_float_set
#define gsl_matrix_TYPE_get gsl_matrix_float_get
#define gsl_matrix_TYPE_set_zero gsl_matrix_float_set_zero
#define gsl_matrix_TYPE_add gsl_matrix_float_add
#define gsl_matrix_TYPE_scale gsl_matrix_float_scale
#define gsl_matrix_TYPE_free gsl_matrix_float_free
#define gsl_matrix_TYPE_view gsl_matrix_float_view
#define gsl_matrix_TYPE_row gsl_matrix_float_row
#define gsl_matrix_TYPE_column gsl_matrix_float_column
#define gsl_matrix_TYPE_get_row gsl_matrix_float_get_row
#define gsl_matrix_TYPE_get_col gsl_matrix_float_get_col
#define gsl_matrix_TYPE_set_row gsl_matrix_float_set_row
#define gsl_matrix_TYPE_set_col gsl_matrix_float_set_col
#define gsl_matrix_TYPE_const_view gsl_matrix_float_const_view
#define gsl_matrix_TYPE_const_submatrix gsl_matrix_float_const_submatrix
#define gsl_matrix_TYPE_fread gsl_matrix_float_fread
#define gsl_matrix_TYPE_fwrite gsl_matrix_float_fwrite
#define gsl_matrix_TYPE_submatrix gsl_matrix_float_submatrix
#define gsl_matrix_TYPE_mul_elements gsl_matrix_float_mul_elements
#define gsl_matrix_TYPE_div_elements gsl_matrix_float_div_elements
#define gsl_matrix_TYPE_add_constant gsl_matrix_float_add_constant
#define gsl_stats_TYPE_min gsl_stats_float_min
#define gsl_stats_TYPE_max gsl_stats_float_max
#define gsl_stats_TYPE_mean gsl_stats_float_mean
#define gsl_stats_TYPE_wmean gsl_stats_float_wmean
#define gsl_stats_TYPE_sd gsl_stats_float_sd
#define gsl_stats_TYPE_median_from_sorted_data gsl_stats_float_median_from_sorted_data
#define gsl_sort_vector_TYPE gsl_sort_vector_float
#define gsl_sort_vector_TYPE_index gsl_sort_vector_float_index
#define gsl_permute_vector_TYPE gsl_permute_vector_float
#define gsl_blas_TYPEcopy gsl_blas_scopy
#define gsl_blas_TYPEaxpy gsl_blas_saxpy
#define gsl_blas_TYPEsyr gsl_blas_ssyr
#define gsl_blas_TYPEtrmv gsl_blas_strmv
#define gsl_blas_TYPEtrmm gsl_blas_strmm
#define gsl_blas_TYPEgemv gsl_blas_sgemv
#define gsl_blas_TYPEgemm gsl_blas_sgemm
#define gsl_blas_TYPEsymm gsl_blas_ssymm
#define gsl_blas_TYPEdot gsl_blas_sdot
#define gsl_blas_TYPEasum gsl_blas_sasum
#define gsl_blas_TYPEnrm2 gsl_blas_snrm2
#define NATIVE_TYPE NATIVE_FLOAT
#endif

#define BLOCKDIM 512

const TYPE const_1=1.0, const_1m=-1.0, const_0=0.0;

#endif
