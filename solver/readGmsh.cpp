#include "readGmsh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>



GmshData readGmsh(const std::string& filepath) {
    std::ifstream file(filepath);
    std::string line;
    std::vector<Node> nodes;
    std::vector<Triangle> triangles;
    std::vector<Quadrangle> quadrangles;

    std::size_t totalNumNodes = 0;

    while (std::getline(file, line)){
        if (line == "$Nodes\r"){
            
            std::getline(file,line); // get first line after nodes
            std::istringstream blk(line);
            std::size_t numEntityBlocks, numNodes, minNodeTag, maxNodeTag;
            blk >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;

            std::cout << "Total number of entity blocks in Gmsh file is: " << numEntityBlocks << "\n";

            for (std::size_t i = 0; i < numEntityBlocks; ++i){ //loop over different entity blocks

                std::getline(file,line);
                std::istringstream blk2(line);
                int entityDim, entityTag, parametric;
                std::size_t numNodesInBlock;
                blk2 >> entityDim >> entityTag >> parametric >> numNodesInBlock ;

                std::cout << "For entity number "<< i+1 << ", the mesh is " << entityDim << "D" << " with " << numNodes << " nodes." << "\n";

                std::vector<std::size_t> nodeTags(numNodesInBlock);
                for (std::size_t j = 0; j<numNodesInBlock; ++j){
                    std::getline(file, line);
                    nodeTags[j] = std::stoull(line);
                }

                for (std::size_t j = 0; j < numNodesInBlock; ++j){
                    std::getline(file, line);
                    std::istringstream coord(line);
                    double x, y, z;
                    coord >> x >> y >> z;

                    nodes.push_back({nodeTags[j], x, y, z});
                }
            
            while (std::getline(file, line)){
                if (line == "$Elements\r"){
                    std::getline(file, line);
                    std::istringstream elem(line);
                    std::size_t numElements, minElementTag, maxElementTag;
                    elem >> numEntityBlocks >> numElements >> minElementTag >> maxElementTag;

                    for (std::size_t i = 0; i< numEntityBlocks; ++i){
                        std::getline(file, line);
                        std::istringstream elem2(line);
                        int elementType;
                        std::size_t numElementsInBlock;

                        elem2 >> entityDim >> entityTag >> elementType >> numElementsInBlock;

                        std::vector<std::size_t> nodeTags(numElementsInBlock);

                        if (elementType == 2){
                            std::cout << "Provided mesh file is 3 node triangle" << "\n";
                            for (std::size_t j = 0; j < numElementsInBlock; ++j){
                                std::getline(file,line);
                                std::istringstream elem3(line);
                                std::size_t x, y, z;
                                elem3 >> x >> y >> z;
    
                                triangles.push_back({j, x, y, z});
                            }
                        }
                        if (elementType == 3){
                            std::cout << "Provided mesh file is 4 node quadrangle" << "\n";
                            for (std::size_t j = 0; j < numElementsInBlock; ++j){
                                std::getline(file,line);
                                std::istringstream elem3(line);
                                std::size_t x, y, z, t;
                                elem3 >> x >> y >> z >> t;
    
                                quadrangles.push_back({j, x, y, z, t});
                            }
                        }
                    }
                }
            }
            }            
        }       
    }
    return GmshData{nodes, triangles, quadrangles};
}
