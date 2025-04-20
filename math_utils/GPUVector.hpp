#ifndef MATH_UTILS_VECTOR_JIM
#define MATH_UTILS_VECTOR_JIM

#include <vector>
#include <Kokkos_Core.hpp>

class GPUVector {
public:
    GPUVector(const std::vector<double>& vector);

    GPUVector(const std::size_t& num_elements);

    std::size_t size() const { return size_; }

    Kokkos::View<double*> get_data_device() const { return data_; }

private:
    Kokkos::View<double*> data_;
    std::size_t size_;
};

#endif
