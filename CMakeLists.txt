cmake_minimum_required(VERSION 3.14.0)
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(fp)

if(NOT DEFINED CMAKE_MODULE_LINKER_FLAGS_MAINTAINER)
    set(
        CMAKE_MODULE_LINKER_FLAGS_MAINTAINER ""
        CACHE STRING ""
    )
endif()

Include(FetchContent)

FetchContent_Declare(
  lemon
  URL http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz
)
FetchContent_MakeAvailable(lemon)

FetchContent_Declare(
    Pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
)
FetchContent_MakeAvailable(Pybind11)

INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_PATH})
include_directories(${lemon_BINARY_DIR})
include_directories(${lemon_SOURCE_DIR})
pybind11_add_module(lemonpy "src/flamingpy/lemonpy/lemonpy.cpp")
target_link_libraries(lemonpy PRIVATE lemon)

add_compile_options(-O3)
