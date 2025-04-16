#ifndef MATH_UTILS_MATRIX_JIM
#define MATH_UTILS_MATRIX_JIM
#include <Kokkos_Core.hpp>
#include "GPUVector.hpp"
#include <vector>
#include <algorithm>
#include <stdexcept>
/**
 * Stores a 2D matrix in CSR format
 * Can be constructed from a 2D matrix std::vector<std::vector<double>> and the values
 * will be inserted. Otherwise an array of {i,j} pairs can be provided which are the locations
 * of non-zero values. They can be updated with update_val(). 
 */
class CSRMatrix {
    public:

        CSRMatrix(const std::vector<std::vector<double>>& matrix);

        CSRMatrix(std::vector<std::pair<int, int>>& nonzero_indices, std::size_t& num_rows, std::size_t& num_cols);

        void update_val(const std::size_t& i, const std::size_t& j, const double& val);

        void copy_to_device();

        GPUVector vmultiply(const GPUVector& x_device) const;

        Kokkos::View<int*> get_rowptr_device(){return row_ptr_;}
        Kokkos::View<int*> get_colind_device(){return col_indices_;}
        Kokkos::View<double*> get_vals_device(){return values_;}
        std::vector<int> get_rowptr_host(){return row_ptr_host;}
        std::vector<int> get_colind_host(){return col_indices_host;}
        std::vector<double> get_vals_host(){return values_host;}

    private:
        Kokkos::View<int*> row_ptr_;
        Kokkos::View<int*> col_indices_;
        Kokkos::View<double*> values_;
        std::vector<int> row_ptr_host;
        std::vector<int> col_indices_host;
        std::vector<double> values_host;
        std::size_t num_rows_ = 0;
        std::size_t num_cols_ = 0;
        std::size_t num_nz_ = 0;
};
#endif