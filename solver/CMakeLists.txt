add_library(solver readGmsh.cpp readKappa.cpp)

# Include current dir for solver headers
target_include_directories(solver PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# Link solver to math_utils and Kokkos
target_link_libraries(solver
    PUBLIC math_utils
    PRIVATE Kokkos::kokkos
)