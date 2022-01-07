# Copyright 2021 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import sys
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

# The following two classes are adaptations of the Python example for pybind11:
# https://github.com/pybind/python_example/blob/master/setup.py

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir="ft_stack/lemonpy"):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        path = self.get_ext_fullpath(ext.name)
        extdir = os.path.abspath(os.path.dirname(path))

        # Required for auto-detection of auxiliary "native" libs.
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        os.makedirs(self.build_temp, exist_ok=True)

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        cmake_args = [
            "-DCMAKE_BUILD_TYPE=Release",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
        ]

        cmake_cmd = ["cmake", ext.sourcedir] + cmake_args
        build_cmd = ["cmake", "--build", "."]

        subprocess.check_call(cmake_cmd, cwd=self.build_temp)
        subprocess.check_call(build_cmd, cwd=self.build_temp)

setup(
    name="ft-stack",
    version="0.1.0",
    description="Threshold estimations for concatenated quantum codes",
    url="https://github.com/XanaduAI/ft-stack",
    packages=find_packages(),
    package_data={"ft_stack":["data/*"]},
    cmdclass={"build_ext": CMakeBuild},
    ext_modules=[CMakeExtension(name="lemonpy")],
    install_requires=[
        "matplotlib==3.3.3",
        "networkx==2.5",
        "retworkx==0.10.2",
        "numpy==1.21.0",
        "pandas==1.2.1",
        "scipy==1.6.0",
        "thewalrus==0.15.0",
        "cmake",
        "pyqt5"
    ]
)
