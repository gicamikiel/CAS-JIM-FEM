add_library(math_utils vector.cpp)

target_include_directories(math_utils PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(math_utils PRIVATE Kokkos::kokkos)