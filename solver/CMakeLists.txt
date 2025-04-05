add_library(solver mesh.cpp)
target_include_directories(solver PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(solver PUBLIC math_utils)