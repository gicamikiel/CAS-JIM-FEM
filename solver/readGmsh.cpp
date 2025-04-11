#include "readGmsh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <Kokkos_Core.hpp>


GmshDataKokkos readGmsh(const std::string& filepath) {
    std::ifstream file(filepath);
    std::string line;

    std::vector<std::size_t> hostNodeTags;
    std::vector<double> hostNodeCoords; // flattened (x, y, z)

    std::vector<std::size_t> hostTriangles;   // flattened (n0, n1, n2)
    std::vector<std::size_t> hostQuadrangles; // flattened (n0, n1, n2, n3)

    while (std::getline(file, line)) {
        if (line == "$Nodes\r") {
            std::getline(file, line); // header
            std::istringstream blk(line);
            std::size_t numEntityBlocks, numNodes, minNodeTag, maxNodeTag;
            blk >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;

            for (std::size_t i = 0; i < numEntityBlocks; ++i) {
                std::getline(file, line);
                std::istringstream blk2(line);
                int entityDim, entityTag, parametric;
                std::size_t numNodesInBlock;
                blk2 >> entityDim >> entityTag >> parametric >> numNodesInBlock;

                std::vector<std::size_t> nodeTags(numNodesInBlock);
                for (std::size_t j = 0; j < numNodesInBlock; ++j) {
                    std::getline(file, line);
                    nodeTags[j] = std::stoull(line);
                }

                for (std::size_t j = 0; j < numNodesInBlock; ++j) {
                    std::getline(file, line);
                    std::istringstream coord(line);
                    double x, y, z;
                    coord >> x >> y >> z;

                    hostNodeTags.push_back(nodeTags[j]);
                    hostNodeCoords.push_back(x);
                    hostNodeCoords.push_back(y);
                    hostNodeCoords.push_back(z);
                }
            }
        }

        if (line == "$Elements\r") {
            std::getline(file, line);
            std::istringstream elem(line);
            std::size_t numEntityBlocks, numElements, minElementTag, maxElementTag;
            elem >> numEntityBlocks >> numElements >> minElementTag >> maxElementTag;

            for (std::size_t i = 0; i < numEntityBlocks; ++i) {
                std::getline(file, line);
                std::istringstream elem2(line);
                int entityDim, entityTag, elementType;
                std::size_t numElementsInBlock;
                elem2 >> entityDim >> entityTag >> elementType >> numElementsInBlock;

                for (std::size_t j = 0; j < numElementsInBlock; ++j) {
                    std::getline(file, line);
                    std::istringstream elem3(line);
                    if (elementType == 2) {
                        // 3-node triangle
                        std::size_t a, b, c;
                        elem3 >> a >> b >> c;
                        hostTriangles.push_back(a);
                        hostTriangles.push_back(b);
                        hostTriangles.push_back(c);
                    } else if (elementType == 3) {
                        // 4-node quad
                        std::size_t a, b, c, d;
                        elem3 >> a >> b >> c >> d;
                        hostQuadrangles.push_back(a);
                        hostQuadrangles.push_back(b);
                        hostQuadrangles.push_back(c);
                        hostQuadrangles.push_back(d);
                    }
                }
            }
        }
    }

    // ===============================
    // Allocate Kokkos Views and copy
    // ===============================

    std::size_t numNodes = hostNodeTags.size();
    std::size_t numTriangles = hostTriangles.size() / 3;
    std::size_t numQuadrangles = hostQuadrangles.size() / 4;

    Kokkos::View<std::size_t*> nodeTags("nodeTags", numNodes);
    Kokkos::View<double**> nodeCoords("nodeCoords", numNodes, 3);

    Kokkos::View<std::size_t*[3]> triangles("triangles", numTriangles);
    Kokkos::View<std::size_t*[4]> quadrangles("quadrangles", numQuadrangles);

    // Host mirrors for copying
    auto host_nodeTags = Kokkos::create_mirror_view(nodeTags);
    auto host_nodeCoords = Kokkos::create_mirror_view(nodeCoords);

    auto host_triangles = Kokkos::create_mirror_view(triangles);
    auto host_quadrangles = Kokkos::create_mirror_view(quadrangles);

    for (std::size_t i = 0; i < numNodes; ++i) {
        host_nodeTags(i) = hostNodeTags[i];
        host_nodeCoords(i, 0) = hostNodeCoords[3 * i + 0];
        host_nodeCoords(i, 1) = hostNodeCoords[3 * i + 1];
        host_nodeCoords(i, 2) = hostNodeCoords[3 * i + 2];
    }

    for (std::size_t i = 0; i < numTriangles; ++i) {
        for (int j = 0; j < 3; ++j)
            host_triangles(i, j) = hostTriangles[3 * i + j];
    }

    for (std::size_t i = 0; i < numQuadrangles; ++i) {
        for (int j = 0; j < 4; ++j)
            host_quadrangles(i, j) = hostQuadrangles[4 * i + j];
    }

    // Deep copy to device
    Kokkos::deep_copy(nodeTags, host_nodeTags);
    Kokkos::deep_copy(nodeCoords, host_nodeCoords);
    Kokkos::deep_copy(triangles, host_triangles);
    Kokkos::deep_copy(quadrangles, host_quadrangles);

    return GmshDataKokkos{
        nodeCoords,
        nodeTags,
        triangles,
        quadrangles
    };
}
