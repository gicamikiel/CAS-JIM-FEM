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
/* 
struct Node {
    std::size_t tag;
    double x, y, z;
};

struct Triangle{
    std::size_t tag;
    std::size_t x, y, z;
};

struct Quadrangle{
    std::size_t tag;
    std::size_t x, y, z, t;
};

struct GmshData{
    std::vector<Node> nodes;
    std::vector<Triangle> triangles;
    std::vector<Quadrangle> quadrangles;
}; */

GmshDataKokkos readGmsh(const std::string& filepath);

#endif