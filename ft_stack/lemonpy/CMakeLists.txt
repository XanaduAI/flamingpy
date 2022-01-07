cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(lemonpy)

if(NOT DEFINED CMAKE_MODULE_LINKER_FLAGS_MAINTAINER)
    set(
        CMAKE_MODULE_LINKER_FLAGS_MAINTAINER ""
        CACHE STRING ""
    )
endif()

Include(FetchContent)

FetchContent_Declare(
    Pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
)

FetchContent_MakeAvailable(Pybind11)

FetchContent_Declare(
  lemon
  URL http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz
)

FetchContent_MakeAvailable(lemon)

include_directories(${lemon_BINARY_DIR})
include_directories(${lemon_SOURCE_DIR})
pybind11_add_module(lemonpy lemonpy.cpp)
target_link_libraries(lemonpy PRIVATE lemon)
