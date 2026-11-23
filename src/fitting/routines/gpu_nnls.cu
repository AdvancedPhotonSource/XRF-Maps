#include "gpu_nnls.hpp"

#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>

#define CUDA_CHECK(call)                                                         \
    do {                                                                         \
        cudaError_t err = (call);                                                \
        if (err != cudaSuccess) {                                                \
            throw std::runtime_error(std::string("CUDA error: ") +               \
                                     cudaGetErrorString(err));                   \
        }                                                                        \
    } while (0)

#define CUBLAS_CHECK(call)                                                       \
    do {                                                                         \
        cublasStatus_t status = (call);                                          \
        if (status != CUBLAS_STATUS_SUCCESS) {                                   \
            throw std::runtime_error("cuBLAS error");                            \
        }                                                                        \
    } while (0)

namespace gpu_nnls {

__global__ void projected_gradient_step_kernel(
    const double* x,
    const double* grad,
    double* x_new,
    double alpha,
    int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        double val = x[i] - alpha * grad[i];
        x_new[i] = (val < 0.0) ? 0.0 : val;
    }
}

__global__ void residual_inf_norm_kernel(
    const double* x,
    const double* grad,
    double* block_max,
    int n)
{
    extern __shared__ double sdata[];
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    double v = 0.0;
    if (i < n) {
        if (x[i] > 0.0) {
            v = fabs(grad[i]);
        } else {
            v = fmax(0.0, -grad[i]);
        }
    }

    sdata[tid] = v;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] = fmax(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }

    if (tid == 0) {
        block_max[blockIdx.x] = sdata[0];
    }
}

GpuNnlsSolver::GpuNnlsSolver() : handle_(nullptr) {
    CUBLAS_CHECK(cublasCreate(&handle_));
}

GpuNnlsSolver::~GpuNnlsSolver() {
    if (handle_ != nullptr) {
        cublasDestroy(handle_);
        handle_ = nullptr;
    }
}

void GpuNnlsSolver::solve(const double* d_A,
                          const double* d_b,
                          int m,
                          int n,
                          double* d_x_out,
                          int max_iters,
                          double tol)
{
    if (m <= 0 || n <= 0) {
        throw std::invalid_argument("Invalid matrix dimensions");
    }
    if (d_A == nullptr || d_b == nullptr || d_x_out == nullptr) {
        throw std::invalid_argument("Null device pointer passed to solve()");
    }

    double *d_G = nullptr, *d_c = nullptr, *d_grad = nullptr;
    double *d_x0 = nullptr, *d_x1 = nullptr;
    double *d_block_max = nullptr;

    try {
        CUDA_CHECK(cudaMalloc(&d_G, static_cast<size_t>(n) * n * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_c, static_cast<size_t>(n) * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_grad, static_cast<size_t>(n) * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_x0, static_cast<size_t>(n) * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_x1, static_cast<size_t>(n) * sizeof(double)));

        const double one = 1.0;
        const double zero = 0.0;
        const double minus_one = -1.0;

        // G = A^T A, upper triangle only.
        CUBLAS_CHECK(cublasDsyrk(
            handle_,
            CUBLAS_FILL_MODE_UPPER,
            CUBLAS_OP_T,
            n,
            m,
            &one,
            d_A,
            m,
            &zero,
            d_G,
            n));

        // Symmetrize on host for simplicity.
        std::vector<double> h_G(static_cast<size_t>(n) * n, 0.0);
        CUDA_CHECK(cudaMemcpy(
            h_G.data(), d_G, static_cast<size_t>(n) * n * sizeof(double), cudaMemcpyDeviceToHost));

        for (int col = 0; col < n; ++col) {
            for (int row = 0; row < n; ++row) {
                if (row > col) {
                    h_G[col * n + row] = h_G[row * n + col];
                }
            }
        }

        CUDA_CHECK(cudaMemcpy(
            d_G, h_G.data(), static_cast<size_t>(n) * n * sizeof(double), cudaMemcpyHostToDevice));

        // c = A^T b
        CUBLAS_CHECK(cublasDgemv(
            handle_,
            CUBLAS_OP_T,
            m,
            n,
            &one,
            d_A,
            m,
            d_b,
            1,
            &zero,
            d_c,
            1));

        CUDA_CHECK(cudaMemset(d_x0, 0, static_cast<size_t>(n) * sizeof(double)));

        // Conservative step size.
        double G_inf = 0.0;
        for (int i = 0; i < n; ++i) {
            double row_sum = 0.0;
            for (int j = 0; j < n; ++j) {
                row_sum += std::abs(h_G[j * n + i]);
            }
            G_inf = std::max(G_inf, row_sum);
        }
        const double alpha = (G_inf > 0.0) ? (1.0 / G_inf) : 1.0;

        const int threads = 256;
        const int blocks = (n + threads - 1) / threads;

        CUDA_CHECK(cudaMalloc(&d_block_max, static_cast<size_t>(blocks) * sizeof(double)));
        std::vector<double> h_block_max(blocks, 0.0);

        double* d_x_curr = d_x0;
        double* d_x_next = d_x1;

        for (int iter = 0; iter < max_iters; ++iter) {
            // grad = G*x_curr - c
            CUBLAS_CHECK(cublasDgemv(
                handle_,
                CUBLAS_OP_N,
                n,
                n,
                &one,
                d_G,
                n,
                d_x_curr,
                1,
                &zero,
                d_grad,
                1));

            CUBLAS_CHECK(cublasDaxpy(
                handle_,
                n,
                &minus_one,
                d_c,
                1,
                d_grad,
                1));

            projected_gradient_step_kernel<<<blocks, threads>>>(
                d_x_curr, d_grad, d_x_next, alpha, n);
            CUDA_CHECK(cudaGetLastError());

            // Compute gradient at x_next for KKT test.
            CUBLAS_CHECK(cublasDgemv(
                handle_,
                CUBLAS_OP_N,
                n,
                n,
                &one,
                d_G,
                n,
                d_x_next,
                1,
                &zero,
                d_grad,
                1));

            CUBLAS_CHECK(cublasDaxpy(
                handle_,
                n,
                &minus_one,
                d_c,
                1,
                d_grad,
                1));

            residual_inf_norm_kernel<<<blocks, threads, threads * sizeof(double)>>>(
                d_x_next, d_grad, d_block_max, n);
            CUDA_CHECK(cudaGetLastError());

            CUDA_CHECK(cudaMemcpy(
                h_block_max.data(),
                d_block_max,
                static_cast<size_t>(blocks) * sizeof(double),
                cudaMemcpyDeviceToHost));

            double kkt_inf = 0.0;
            for (double v : h_block_max) {
                kkt_inf = std::max(kkt_inf, v);
            }

            if (kkt_inf < tol) {
                CUDA_CHECK(cudaMemcpy(
                    d_x_out, d_x_next, static_cast<size_t>(n) * sizeof(double), cudaMemcpyDeviceToDevice));
                std::cout << "NNLS converged in " << (iter + 1)
                          << " iterations, KKT residual = " << kkt_inf << "\n";
                break;
            }

            std::swap(d_x_curr, d_x_next);

            if (iter == max_iters - 1) {
                CUDA_CHECK(cudaMemcpy(
                    d_x_out, d_x_curr, static_cast<size_t>(n) * sizeof(double), cudaMemcpyDeviceToDevice));
                std::cout << "NNLS reached max iterations, KKT residual = "
                          << kkt_inf << "\n";
            }
        }

        CUDA_CHECK(cudaFree(d_G));
        CUDA_CHECK(cudaFree(d_c));
        CUDA_CHECK(cudaFree(d_grad));
        CUDA_CHECK(cudaFree(d_x0));
        CUDA_CHECK(cudaFree(d_x1));
        CUDA_CHECK(cudaFree(d_block_max));
    } catch (...) {
        if (d_G) cudaFree(d_G);
        if (d_c) cudaFree(d_c);
        if (d_grad) cudaFree(d_grad);
        if (d_x0) cudaFree(d_x0);
        if (d_x1) cudaFree(d_x1);
        if (d_block_max) cudaFree(d_block_max);
        throw;
    }
}

} // namespace gpu_nnls

