cmake_minimum_required(VERSION 3.23)
project(GMISolver)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Werror -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Werror")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(GUROBI REQUIRED)
include(get_cpm.cmake)

CPMAddPackage("gh:gabime/spdlog@1.10.0")
CPMAddPackage("gh:fmtlib/fmt#9.1.0")
CPMAddPackage("gh:abseil/abseil-cpp#20220623.1")

set(BASEPATH "${CMAKE_SOURCE_DIR}")
include_directories("${BASEPATH}")
include_directories("${GUROBI_INCLUDE_DIRS}")
add_subdirectory(src)
add_subdirectory(tests)

add_executable(GMISolver main.cpp)
target_link_libraries(GMISolver
        Algorithms
        DataModel
        Util
        fmt::fmt-header-only
        spdlog::spdlog
        absl::flags
        absl::flags_parse)
target_link_libraries(GMISolver optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY})
target_link_libraries(GMISolver ${GUROBI_LIBRARY})
