add_executable(unit_tests linalg.cpp)
target_link_libraries(unit_tests PRIVATE Catch2::Catch2WithMain)
target_link_libraries(unit_tests PUBLIC solver)

include(CTest)
include(Catch)

catch_discover_tests(unit_tests)