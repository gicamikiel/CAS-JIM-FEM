add_executable(unit_tests unitTests.cpp linalg.cpp reader.cpp shapeFunctions.cpp)

target_link_libraries(unit_tests PRIVATE Catch2::Catch2)
target_link_libraries(unit_tests PRIVATE solver)
target_link_libraries(unit_tests PRIVATE Kokkos::kokkos)

include(CTest)
include(Catch)

catch_discover_tests(unit_tests)