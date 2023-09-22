include(../get_cpm.cmake)
cpmaddpackage(
        NAME
        googletest
        GITHUB_REPOSITORY
        google/googletest
        GIT_TAG
        release-1.12.1
        VERSION
        1.12.1
        OPTIONS
        "INSTALL_GTEST OFF"
        "gtest_force_shared_crt")

add_library(CommonDefs CommonDefs.cpp CommonDefs.h)
target_link_libraries(CommonDefs EnumTypes absl::flags)

add_executable(UnitTests DualSimplexTest.cpp PrimalSimplexTest.cpp)
include_directories(../src/)
target_link_libraries(
        UnitTests
        gtest
        gtest_main
        gmock
        CommonDefs
        RevisedDualSimplexPFIBounds
        RevisedPrimalSimplexPFIBounds
        SimplexTableau
        GurobiOptimizer
        LPOptStatistics
        MpsReader)
enable_testing()
add_test(NAME UnitTests COMMAND UnitTests)
