#include <catch2/catch_test_macros.hpp>
#include <Kokkos_Core.hpp>
#include "GPUVector.hpp"
#include "CSRMatrix.hpp"
#include <vector>
#include <iostream>

TEST_CASE("GPUVector creation", "[GPUVector]") {
    std::vector<double> hostVec = {1.0, 2.0, 3.0, 4.0};
    GPUVector gpuVec(hostVec);

    auto mirror = Kokkos::create_mirror_view(gpuVec.get_data_device());
    Kokkos::deep_copy(mirror, gpuVec.get_data_device());

    for (std::size_t i = 0; i < hostVec.size(); ++i) {
        REQUIRE(mirror(i) == hostVec[i]);
    }
}


TEST_CASE("CSRMatrix creation from sparse matrix", "[CSRMatrix]") {
    // Create a sample sparse matrix.
    // Matrix layout:
    // [ 0   0   0   0   1]
    // [ 5   0   0   8   0]
    // [ 0   0   3   0   0]
    // [ 0   6   0   0   1]
    std::vector<std::vector<double>> sparseMatrix = {
        {0.0, 0.0, 0.0, 0.0, 1.0},
        {5.0, 0.0, 0.0, 8.0, 0.0},
        {0.0, 0.0, 3.0, 0.0, 0.0},
        {0.0, 6.0, 0.0, 0.0, 1.0}
    };

    CSRMatrix csr(sparseMatrix);
    csr.copy_to_device();

    std::vector<int> expectedRowPtr     = {0, 1, 3, 4, 6};
    std::vector<int> expectedColIndices = {4, 0, 3, 2, 1, 4};
    std::vector<double> expectedValues  = {1.0, 5.0, 8.0, 3.0, 6.0, 1.0};

    // Compare CSRMatrix's host arrays with the expected values.
    auto mirror_rowptr = Kokkos::create_mirror_view(csr.get_rowptr_device());
    auto mirror_colind = Kokkos::create_mirror_view(csr.get_colind_device());
    auto mirror_vals = Kokkos::create_mirror_view(csr.get_vals_device());

    REQUIRE(csr.get_rowptr_host() == expectedRowPtr);
    REQUIRE(csr.get_colind_host() == expectedColIndices);
    REQUIRE(csr.get_vals_host() == expectedValues);
}

TEST_CASE("CSRMatrix creation from nonzero indices", "[CSRMatrix]") {
    // Create a sample sparse matrix.
    // Matrix layout:
    // [ 0   0   0   0   1]
    // [ 5   0   0   8   0]
    // [ 0   0   3   0   0]
    // [ 0   6   0   0   1]
    std::vector<std::pair<int, int>> nz_indices = {
        {0,4},
        {1,0}, {1,3},
        {2,2},
        {3,1}, {3,4}
    };

    CSRMatrix csr(nz_indices,4,5);
    csr.copy_to_device();

    std::vector<int> expectedRowPtr     = {0, 1, 3, 4, 6};
    std::vector<int> expectedColIndices = {4, 0, 3, 2, 1, 4};
    std::vector<double> expectedValues  = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    // Compare CSRMatrix's host arrays with the expected values.
    auto mirror_rowptr = Kokkos::create_mirror_view(csr.get_rowptr_device());
    auto mirror_colind = Kokkos::create_mirror_view(csr.get_colind_device());
    auto mirror_vals = Kokkos::create_mirror_view(csr.get_vals_device());

    REQUIRE(csr.get_rowptr_host() == expectedRowPtr);
    REQUIRE(csr.get_colind_host() == expectedColIndices);
    REQUIRE(csr.get_vals_host() == expectedValues);
}

TEST_CASE("CSRMatrix update values on host before copy", "[CSRMatrix]") {
    // Create a sample sparse matrix.
    // Matrix layout:
    // [ 0   0   0   0   1]
    // [ 5   0   0   8   0]
    // [ 0   0   3   0   0]
    // [ 0   6   0   0   1]

    //Start with "empty" matrix
    std::vector<std::pair<int, int>> nz_indices = {
        {0,4},
        {1,0}, {1,3},
        {2,2},
        {3,1}, {3,4}
    };

    CSRMatrix csr(nz_indices,4,5);
    csr.update_val(2,2,3.0);
    csr.copy_to_device();

    std::vector<int> expectedRowPtr     = {0, 1, 3, 4, 6};
    std::vector<int> expectedColIndices = {4, 0, 3, 2, 1, 4};
    std::vector<double> expectedValues  = {1.0, 1.0, 1.0, 3.0, 1.0, 1.0};

    // Compare CSRMatrix's host arrays with the expected values.
    auto mirror_rowptr = Kokkos::create_mirror_view(csr.get_rowptr_device());
    auto mirror_colind = Kokkos::create_mirror_view(csr.get_colind_device());
    auto mirror_vals = Kokkos::create_mirror_view(csr.get_vals_device());

    REQUIRE(csr.get_rowptr_host() == expectedRowPtr);
    REQUIRE(csr.get_colind_host() == expectedColIndices);
    REQUIRE(csr.get_vals_host() == expectedValues);
}


TEST_CASE("CSRMatrix and GPUVector multiplication", "[CSRMatrix][GPUVector]") {
    std::vector<std::vector<double>> sparseMatrix = {
        {0.0, 0.0, 0.0, 0.0, 1.0},
        {5.0, 0.0, 0.0, 8.0, 0.0},
        {0.0, 0.0, 3.0, 0.0, 0.0},
        {0.0, 6.0, 0.0, 0.0, 1.0}
    };

    CSRMatrix csr(sparseMatrix);
    csr.copy_to_device();

    std::vector<double> hostVec = {1.0, 2.0, 3.0, 4.0, 5.0};
    GPUVector gpuVec(hostVec);

    // Perform GPU matrix-vector multiplication
    GPUVector y_device = csr.vmultiply(gpuVec);

    auto y_view = y_device.get_data_device();
    std::cout <<"Breakpoint 2" << '\n';
    auto y_host = Kokkos::create_mirror_view(y_view);
    Kokkos::deep_copy(y_host, y_view);

    std::vector<double> expected = {5.0, 37.0, 9.0, 17.0};

    for (std::size_t i = 0; i < expected.size(); ++i) {
        CAPTURE(i);
        REQUIRE(y_host(i) == expected[i]);
    }
}
