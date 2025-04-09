#ifndef READ_GMSH
#define READ_GMSH
#include <vector>
#include <string>

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
};

GmshData readGmsh(const std::string& filepath);

#endif