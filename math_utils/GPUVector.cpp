#include <GPUVector.hpp>

GPUVector::GPUVector(const std::vector<double>& vector) {
    size_ = vector.size();
    data_ = Kokkos::View<double*>("data", size_);
    
    auto host_mirror = Kokkos::create_mirror_view(data_);

    for (std::size_t i = 0; i < size_; ++i) {
        host_mirror(i) = vector[i];
    }
    Kokkos::deep_copy(data_, host_mirror);
}

GPUVector::GPUVector(const std::size_t& num_elements) {
    size_ = num_elements;
    data_ = Kokkos::View<double*>("data", size_);
    
    auto host_mirror = Kokkos::create_mirror_view(data_);
    Kokkos::deep_copy(data_, host_mirror);
}



    
