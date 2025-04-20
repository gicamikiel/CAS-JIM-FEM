#include "readKappa.hpp"
#include <Kokkos_Core.hpp>
#include <fstream>
#include <iostream>
#include <sstream>  // For std::istringstream
#include <vector>

KappaDataKokkos readKappa(const std::string& filepath) {
    std::ifstream file(filepath);
    std::string line;

    std::vector<std::size_t> hostElemTags;
    std::vector<double> hostKappas;

    while (std::getline(file, line)) {
        std::istringstream elem(line);
        std::size_t tag; 
        double kappa;
        elem >> tag >> kappa;
        hostElemTags.push_back(tag);
        hostKappas.push_back(kappa);
    }

    std::size_t N = hostElemTags.size();

    Kokkos::View<std::size_t*> elemTags("elemTags", N);
    Kokkos::View<double*> kappas("kappas", N);

    auto h_elemTags = Kokkos::create_mirror_view(elemTags);
    auto h_kappas = Kokkos::create_mirror_view(kappas);

    for (std::size_t i = 0; i < N; ++i) {
        h_elemTags(i) = hostElemTags[i];
        h_kappas(i) = hostKappas[i];
    }

    Kokkos::deep_copy(elemTags, h_elemTags);
    Kokkos::deep_copy(kappas, h_kappas);

    return KappaDataKokkos{elemTags, kappas};
}
