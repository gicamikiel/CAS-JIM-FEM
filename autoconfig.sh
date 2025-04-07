BUILD_TYPE=${1:-Release}
CXX_COMPILER=${2:-g++}
# Backend choices, defaults are for CCI
# https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html#keywords-backends
HOST_BACKEND=${3:THREADS}
DEVICE_BACKEND=${4:CUDA}

if [ ! -d ../build-Catch2 ]; then
    git clone https://github.com/catchorg/Catch2 ../Catch2
    cmake -S ../Catch2 -B ../build/Catch2 \
        -DBUILD_TESTING=OFF \
        -DCMAKE_INSTALL_PREFIX=../build-Catch2/install/
    cmake --build ../build/Catch2 -j$(nproc) --target install
fi

if [ ! -d ../build-Kokkos ]; then
    cd ..
    export KOKKOS_VERSION=4.6.00
    export KOKKOS_DOWNLOAD_URL=https://github.com/kokkos/kokkos/releases/download/${KOKKOS_VERSION}
    curl -sLO ${KOKKOS_DOWNLOAD_URL}/kokkos-${KOKKOS_VERSION}.tar.gz
    curl -sLO ${KOKKOS_DOWNLOAD_URL}/kokkos-${KOKKOS_VERSION}-SHA-256.txt
    grep kokkos-${KOKKOS_VERSION}.tar.gz kokkos-${KOKKOS_VERSION}-SHA-256.txt | shasum -c
    tar -xzvf kokkos-${KOKKOS_VERSION}.tar.gz
    cd -
    cmake -S ../kokkos-${KOKKOS_VERSION} -B ../build/Kokkos \
        -DBUILD_TESTING=OFF \
        -DKokkos_ENABLE_$HOST_BACKEND=ON \
        -DKokkos_ENABLE_$DEVICE_BACKEND=ON \
        -DCMAKE_INSTALL_PREFIX=../build-Kokkos/install/
    cmake --build ../build/Kokkos -j$(nproc) --target install
fi

export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:../build-Catch2/install
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:../build-Kokkos/install
cmake -S . -B ./build \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DKokkos_ENABLE_$HOST_BACKEND=ON \
    -DKokkos_ENABLE_$DEVICE_BACKEND=ON \
    -DCMAKE_CXX_COMPILER=$CXX_COMPILER
cmake --build ./build -j$(nproc)