#ifndef READ_GMSH
#define READ_GMSH
#include <vector>
#include <string>
#include <Kokkos_Core.hpp>

struct GmshDataKokkos {
    Kokkos::View<double**> nodeCoords;
    Kokkos::View<std::size_t*> nodeTags;
    Kokkos::View<std::size_t*[3]> triangles;
    Kokkos::View<std::size_t*[4]> quadrangles;
};

GmshDataKokkos readGmsh(const std::string& filepath);

#endif