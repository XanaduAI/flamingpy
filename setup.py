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
"""setup.py instructions for FlamingPy installation from Source
"""

# pylint: disable=too-few-public-methods,B607

import os
import re
import sys
import platform
import subprocess

from distutils.version import LooseVersion

from setuptools import setup, Extension, dist, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.build_ext import build_ext


# Reading the package version number
with open("flamingpy/_version.py", encoding="utf8") as f:
    version = f.readlines()[-1].split()[-1].strip("\"'")


class BinaryDistribution(dist.Distribution):
    """A class to define Binary Distribution objects"""

    def has_ext_modules(foo):
        """Check for external modules."""
        return True


class CMakeExtension(Extension):
    """A class to define CMake Extensions.

    Adapted from Python examples for pybind11:
    https://github.com/pybind/python_example/blob/master/setup.py
    """

    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    """A class to define, configure, and test build extensions.

    Adapted from the pymatching package:
    https://github.com/oscarhiggott/PyMatching/blob/master/setup.py
    """

    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError as exc:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            ) from exc

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r"version\s*([\d.]+)", out.decode()).group(1))
            if cmake_version < "3.14.0":
                raise RuntimeError("CMake >= 3.14.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)


# use README.md as long_description
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf8") as f:
    long_description = f.read()


# setup parameters
if sys.argv[1] == "build_cython":
    ext_modules = [
        Extension(
            "flamingpy.cpp.cpp_mc_loop",
            sources=["flamingpy/cpp/cpp_mc_loop.pyx"],
            extra_compile_args=["-O3", "-w"],
            language="c++",
        )
    ]
elif sys.argv[1] == "build_cmake":
    ext_modules = [CMakeExtension("flamingpy.cpp.lemonpy")]
elif sys.argv[1] == "install" or sys.argv[1] == "develop" or sys.argv[1] == "bdist_wheel":
    ext_modules = []
else:
    raise NotImplementedError

classifiers = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: Apache Software License",
    "Natural Language :: English",
    "Operating System :: POSIX",
    "Operating System :: MacOS :: MacOS X",
    "Operating System :: POSIX :: Linux",
    "Operating System :: Microsoft :: Windows",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3 :: Only",
    "Topic :: Scientific/Engineering :: Physics",
]

install_requires = [
    "matplotlib>=3.3.3",
    "networkx>=2.5",
    "numpy>=1.21",
    "retworkx>=0.10.2",
    "pandas>=1.2.1",
    "scipy>=1.6",
]

description = """FlamingPy is a cross-platform Python library with a variety of backends for
efficient simulations of error correction in fault-tolerant quantum computers."""

setup(
    name="flamingpy",
    version=version,
    description=description,
    license="Apache License 2.0",
    classifiers=classifiers,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/XanaduAI/flamingpy",
    packages=find_packages(where="."),
    include_package_data=True,
    package_dir={"": "."},
    python_requires=">=3.8,!=3.11.*",
    cmdclass={
        "install": install,
        "develop": develop,
        "build_cython": build_ext,
        "build_cmake": CMakeBuild,
    },
    ext_modules=ext_modules,
    distclass=BinaryDistribution,
    install_requires=install_requires,
)
