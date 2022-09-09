cmake_minimum_required(VERSION 3.23)
project(GMISolver)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Werror -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Werror")

include(get_cpm.cmake)

CPMAddPackage("gh:gabime/spdlog@1.10.0")
CPMAddPackage("gh:fmtlib/fmt#9.1.0")
CPMAddPackage("gh:abseil/abseil-cpp#20220623.1")

set(BASEPATH "${CMAKE_SOURCE_DIR}")
include_directories("${BASEPATH}")
add_subdirectory(src)

add_executable(GMISolver main.cpp)
target_link_libraries(GMISolver
        Algorithms
        DataModel
        Util
        fmt::fmt-header-only
        spdlog::spdlog
        absl::flags
        absl::flags_parse)
