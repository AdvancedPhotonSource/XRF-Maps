#pragma once

#include <cublas_v2.h>

namespace gpu_nnls {

class GpuNnlsSolver {
public:
    GpuNnlsSolver();
    ~GpuNnlsSolver();

    GpuNnlsSolver(const GpuNnlsSolver&) = delete;
    GpuNnlsSolver& operator=(const GpuNnlsSolver&) = delete;

    void solve(const double* d_A,
               const double* d_b,
               int m,
               int n,
               double* d_x_out,
               int max_iters = 5000,
               double tol = 1e-8);

private:
    cublasHandle_t handle_;
};

} // namespace gpu_nnls
