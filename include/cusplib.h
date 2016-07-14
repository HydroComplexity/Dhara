#ifndef __CUSPLIB_H__
#define __CUSPLIB_H__
#endif

#include <cuda_runtime.h>
#include <cusp/array1d.h>
#include <cusp/coo_matrix.h>
#include <cusp/hyb_matrix.h>
#include <cusp/monitor.h>
#include <cusp/multiply.h>
#include <cusp/elementwise.h>
#include <cusp/gallery/poisson.h>
#include <cusp/linear_operator.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/print.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>



///////////////////////////////////////////
// User type definitions                 //
///////////////////////////////////////////

typedef cusp::array1d<double, cusp::device_memory> cuspdev_1d;
typedef cusp::dia_matrix<int, double, cusp::device_memory> cuspdev_diamat;
typedef cusp::identity_operator<double, cusp::device_memory> cuspdev_idoper;
typedef thrust::device_vector<double>::iterator thrustdev_iter;
typedef thrust::device_vector<double> thrustdev;
