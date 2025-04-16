#include "CSRMatrix.hpp"

// Constructor from std::vector<std::vector<double>>
CSRMatrix::CSRMatrix(const std::vector<std::vector<double>>& matrix) {
    num_rows_ = matrix.size();
    num_cols_ = matrix[0].size();
    num_nz_ = 0;
    row_ptr_host.resize(num_rows_ + 1, 0);

    for (std::size_t i = 0; i < num_rows_; ++i) {
        for (std::size_t j = 0; j < num_cols_; ++j) {
            if (matrix[i][j] != 0.0) {
                col_indices_host.push_back(j);
                values_host.push_back(matrix[i][j]);
                ++num_nz_;
            }
        }
        row_ptr_host[i + 1] = num_nz_; 
    }
}

//Constructor from array of i,j non-zero value positions
//values can be changed later with update_val
CSRMatrix::CSRMatrix(std::vector<std::pair<int, int>>& nonzero_indices, std::size_t& num_rows, std::size_t& num_cols) 
: num_rows_(num_rows), num_cols_(num_cols) {
    num_nz_ = nonzero_indices.size();
    row_ptr_host.resize(num_rows_ + 1, 0);
    col_indices_host.resize(num_nz_);
    values_host.resize(num_nz_,1.0);

    // Sort by row-major order to prepare for CSR structure
    std::sort(nonzero_indices.begin(), nonzero_indices.end());

    //TODO: Implement the filling of the sparsity pattern
}

void CSRMatrix::copy_to_device(){
    row_ptr_ = Kokkos::View<int*>("row_ptr", row_ptr_host.size());
    col_indices_ = Kokkos::View<int*>("col_indices", col_indices_host.size());
    values_ = Kokkos::View<double*>("values", values_host.size());

    Kokkos::View<int*, Kokkos::HostSpace> row_ptr_host_view(row_ptr_host.data(), row_ptr_host.size());
    Kokkos::View<int*, Kokkos::HostSpace> col_indices_host_view(col_indices_host.data(), col_indices_host.size());
    Kokkos::View<double*, Kokkos::HostSpace> values_host_view(values_host.data(), values_host.size());

    Kokkos::deep_copy(row_ptr_, row_ptr_host_view);
    Kokkos::deep_copy(col_indices_, col_indices_host_view);
    Kokkos::deep_copy(values_, values_host_view);
}


// Update a value at a specific i, j if it's a valid non-zero
void CSRMatrix::update_val(const std::size_t& i, const std::size_t& j, const double& val) {
    auto start = row_ptr_(i);
    auto end   = row_ptr_(i + 1);
    bool found = false;
    for (int idx = start; idx < end; ++idx) {
        if (col_indices_(idx) == j) {
            values_(idx) = val;
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Specified i,j pair does not exist");
    }
}


GPUVector CSRMatrix::vmultiply(const GPUVector& x_device) const {
    if (x_device.size() != num_cols_) {
        throw std::runtime_error("Matrix and vector dimensions not compatible for multiplication.");
    }

    GPUVector y_device(num_rows_);
    auto x_view = x_device.get_data_device();
    auto y_view = y_device.get_data_device();

    auto row_ptr  = row_ptr_;
    auto col_inds = col_indices_;
    auto values   = values_;

    Kokkos::parallel_for("spmv", num_rows_, KOKKOS_LAMBDA(const std::size_t i) {
        double sum = 0.0;
        for (int idx = row_ptr(i); idx < row_ptr(i + 1); ++idx) {
            sum += values(idx) * x_view(col_inds(idx));
        }
        y_view(i) = sum;
    });

    return y_device;
}








