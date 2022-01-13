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
import re
import sys
import platform
import subprocess

from setuptools import setup, dist, Extension, find_packages
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


# The following class is an adaptation of Python examples for pybind11:
# https://github.com/pybind/python_example/blob/master/setup.py
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

        
# The following class was adapted from pymatching package:
# https://github.com/oscarhiggott/PyMatching/blob/master/setup.py
class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)


class BinaryDistribution(dist.Distribution):
        def has_ext_modules(foo):
                    return True


# use README.md as long_description
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md')) as f:
        long_description = f.read()


setup(
    name="ft-stack",
    version="0.1.12",
    description="Threshold estimations for concatenated quantum codes",
    license="Apache License 2.0",
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/XanaduAI/ft-stack",
    packages=find_packages("src"),
    package_dir={'':'src'},
    #package_data={"ft_stack":["src/ft_stack/data/*.csv", "src/ft_stack/*.so", "examples/lemon_benchmark.py"]},
    include_package_data=True,
    python_requires='>=3.8',
    cmdclass={"build_ext": CMakeBuild},
    ext_modules=[CMakeExtension('ft_stack.lemonpy')],
    distclass=BinaryDistribution,
    install_requires=[
        "matplotlib>=3.3.3",
        "networkx>=2.5",
        "retworkx>=0.10.2",
        "numpy>=1.21.0",
        "pandas>=1.2.1",
        "scipy>=1.6.0",
        "thewalrus>=0.15.0"
    ]
)
