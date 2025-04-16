BUILD_TYPE=${1:-Release}
CXX_COMPILER=${2:-g++}
# Backend choices, defaults are for CCI
# https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html#keywords-backends
HOST_BACKEND=${3:-THREADS}
DEVICE_BACKEND=${4:-CUDA}

if [ ! -d ../Catch2-install ]; then
    git clone https://github.com/catchorg/Catch2 ../Catch2
    rm -rf ../build/Catch2
    cmake -S ../Catch2 -B ../build/Catch2 \
        -DBUILD_TESTING=OFF \
        -DCMAKE_INSTALL_PREFIX=../Catch2-install/
    cmake --build ../build/Catch2 -j$(nproc) --target install
fi

if [ ! -d ../Kokkos-install ]; then
    cd ..
    export KOKKOS_VERSION=4.6.00
    export KOKKOS_DOWNLOAD_URL=https://github.com/kokkos/kokkos/releases/download/${KOKKOS_VERSION}
    curl -sLO ${KOKKOS_DOWNLOAD_URL}/kokkos-${KOKKOS_VERSION}.tar.gz
    curl -sLO ${KOKKOS_DOWNLOAD_URL}/kokkos-${KOKKOS_VERSION}-SHA-256.txt
    grep kokkos-${KOKKOS_VERSION}.tar.gz kokkos-${KOKKOS_VERSION}-SHA-256.txt | shasum -c
    tar -xzvf kokkos-${KOKKOS_VERSION}.tar.gz
    cd -
    rm -rf ../build/Kokkos
    cmake -S ../kokkos-${KOKKOS_VERSION} -B ../build/Kokkos \
        -DBUILD_TESTING=OFF \
        -DKokkos_ENABLE_$HOST_BACKEND=ON \
        -DKokkos_ENABLE_$DEVICE_BACKEND=ON \
        -DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON \
        -DCMAKE_INSTALL_PREFIX=../Kokkos-install/
    cmake --build ../build/Kokkos -j$(nproc) --target install
fi

cmake -S . -B ./build \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DCatch2_ROOT=../Catch2-install \
    -DKokkos_ROOT=../Kokkos-install \
    -DKokkos_ENABLE_$HOST_BACKEND=ON \
    -DKokkos_ENABLE_$DEVICE_BACKEND=ON \
    -DCMAKE_CXX_COMPILER=$CXX_COMPILER
cmake --build ./build -j$(nproc)