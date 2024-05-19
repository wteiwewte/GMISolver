cmake_minimum_required(VERSION 3.23)
project(GMISolver)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_DEBUG
    "-Wall -Werror -Wno-deprecated -Wno-error=array-parameter -fsanitize=address"
)
set(CMAKE_CXX_FLAGS_RELEASE
    "-O3 -DNDEBUG -Wall -Werror -Wno-deprecated -Wno-error=array-parameter")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(GUROBI REQUIRED)
include(get_cpm.cmake)

cpmaddpackage("gh:gabime/spdlog@1.10.0")
cpmaddpackage("gh:fmtlib/fmt#9.1.0")
cpmaddpackage("gh:abseil/abseil-cpp#20220623.1")
cpmaddpackage("gh:TartanLlama/expected@1.1.0")
cpmaddpackage(
  NAME
  Boost
  VERSION
  1.83.0
  GITHUB_REPOSITORY
  "boostorg/boost"
  GIT_TAG
  "boost-1.83.0")

set(BASEPATH "${CMAKE_SOURCE_DIR}")
include_directories("${BASEPATH}")
add_subdirectory(src/Algorithms)
add_subdirectory(src/DataModel)
add_subdirectory(src/Util)
add_subdirectory(tests)

add_executable(GMISolver main.cpp)
target_link_libraries(
  GMISolver
  DualSimplexGomory
  DualSimplex
  PrimalSimplex
  SimplexTableau
  GurobiOptimizer
  LPOptStatisticsPrinter
  MpsReader
  SpdlogHeader
  absl::flags
  absl::flags_parse)
