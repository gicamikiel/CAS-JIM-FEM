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

        //Constructor from array of i,j non-zero value positions
        template<typename IndexT,typename RowCountT,typename ColCountT,typename = std::enable_if_t<std::is_integral_v<IndexT>>>
        CSRMatrix(const std::vector<std::pair<IndexT, IndexT>>& nonzero_indices,RowCountT num_rows,ColCountT num_cols) 
        : num_rows_(static_cast<std::size_t>(num_rows)), num_cols_(static_cast<std::size_t>(num_cols))
        {
            num_nz_ = nonzero_indices.size();

            col_indices_host.resize(num_nz_);
            values_host.resize(num_nz_, 1.0);

            std::vector<std::pair<IndexT,IndexT>> idxs = nonzero_indices;
            std::sort(idxs.begin(), idxs.end());

            std::vector<std::size_t> row_counts(num_rows_, 0);
            for (std::size_t i = 0; i < num_nz_; ++i) {
                auto [r, c] = idxs[i];
                row_counts[static_cast<std::size_t>(r)]++;
                col_indices_host[i] = static_cast<std::size_t>(c);
            }
            row_ptr_host.resize(num_rows_ + 1);
            row_ptr_host[0] = 0;
            for (std::size_t r = 0; r < num_rows_; ++r) {
                row_ptr_host[r+1] = row_ptr_host[r] + row_counts[r];
            }
        }

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