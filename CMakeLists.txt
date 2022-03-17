# Copyright 2022 Xanadu Quantum Technologies Inc.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
cmake_minimum_required(VERSION 3.14)

# Set the project name
project(fp)


# Find required packages
#find_package(PythonLibs REQUIRED)


set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


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


include_directories(SYSTEM ${PYTHON_INCLUDE_PATH})

include_directories(${lemon_BINARY_DIR})
include_directories(${lemon_SOURCE_DIR})
pybind11_add_module(lemonpy "flamingpy/cpp/lemonpy.cpp")
target_link_libraries(lemonpy PRIVATE lemon)
#set_target_properties(lemonpy PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lemonpy)


add_compile_options(-O3 -w)
