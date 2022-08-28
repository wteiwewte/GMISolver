cmake_minimum_required(VERSION 3.23)
project(GMISolver)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Werror -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Werror")

find_package(fmt CONFIG REQUIRED)
Include(FetchContent)
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.9.2
)
FetchContent_MakeAvailable(spdlog)

#FetchContent_Declare(
#        abseil
#        GIT_REPOSITORY https://github.com/abseil/abseil-cpp
#        GIT_TAG 20220623.0
#)
#FetchContent_MakeAvailable(abseil)

set(BASEPATH "${CMAKE_SOURCE_DIR}")
include_directories("${BASEPATH}")
add_subdirectory(src)

add_executable(GMISolver main.cpp)
target_link_libraries(GMISolver Algorithms DataModel Util fmt::fmt-header-only spdlog::spdlog)
