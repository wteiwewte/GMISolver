include(../get_cpm.cmake)
CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        GIT_TAG release-1.12.1
        VERSION 1.12.1
        OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt"
)
CPMAddPackage("gh:cpm-cmake/testpack-fibonacci@2.0")

add_executable(AllTests DualSimplexTest.cpp
        PrimalSimplexTest.cpp
        CommonDefs.cpp
        CommonDefs.h)
include_directories(../src/)
target_link_libraries(AllTests fibonacci gtest gtest_main gmock absl::flags spdlog::spdlog Algorithms Util ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
enable_testing()
add_test(NAME AllTests COMMAND AllTests)