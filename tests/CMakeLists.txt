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
target_link_libraries(CommonDefs EnumTypes absl::flags gtest)

add_library(LPTestBase INTERFACE LPTestBase.h)
target_link_libraries(
  LPTestBase
  INTERFACE CommonDefs
            EnumTypes
            LPOptStatistics
            IPOptStatistics
            OptStatisticsPrinter
            SimplexTableau
            gtest)

add_executable(
  UnitTests DualSimplexTest.cpp PrimalSimplexTest.cpp DualSimplexGomoryTest.cpp
            LexicographicOptimizerTest.cpp)
include_directories(../src/)
target_link_libraries(
  UnitTests
  gtest
  gtest_main
  gmock
  CommonDefs
  DualSimplexGomory
  RevisedDualSimplexPFIBounds
  RevisedPrimalSimplexPFIBounds
  SimplexTableau
  GurobiOptimizer
  LPOptStatistics
  MpsReader
  LPTestBase
  IPOptStatistics)
enable_testing()
add_test(NAME UnitTests COMMAND UnitTests)
