#include <iostream>
//#include "vector.hpp"
#include "readGmsh.hpp"


int main() {
    std::string path = "../../dataset/bun_zipper_res4";
    GmshDataKokkos data = readGmsh(path);

    return 0;
}