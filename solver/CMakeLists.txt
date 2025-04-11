add_library(solver readGmsh.cpp)

# Include current dir for solver headers
target_include_directories(solver PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# Link solver to math_utils and Kokkos
target_link_libraries(solver
    PUBLIC math_utils
    PRIVATE Kokkos::kokkos
)

# Create main executable
add_executable(main main.cpp)

# Properly link main to solver and Kokkos
target_link_libraries(main
    PRIVATE solver
    PRIVATE Kokkos::kokkos
)

# Optional: if main.cpp needs to see solver headers
target_include_directories(main PRIVATE ${CMAKE_CURRENT_LIST_DIR})
