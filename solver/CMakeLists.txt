add_library(solver mesh.cpp readGmsh.cpp)
target_include_directories(solver PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(solver PUBLIC math_utils)
target_link_libraries(solver PRIVATE Kokkos::kokkos)

add_executable(main main.cpp)
target_link_libraries(main PUBLIC solver)
target_include_directories(main PUBLIC ${CMAKE_CURRENT_LIST_DIR})