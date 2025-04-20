#ifndef READ_KAPPA
#define READ_KAPPA
#include <vector>
#include <string>
#include <Kokkos_Core.hpp>

struct KappaDataKokkos {
    Kokkos::View<std::size_t*> elemTags;
    Kokkos::View<double*> kappas;
};

KappaDataKokkos readKappa(const std::string& filepath);

#endif