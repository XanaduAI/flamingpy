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

from setuptools import find_packages, setup

setup(
    name="ft-stack",
    version="0.1.0",
    description="Threshold estimations for concatenated quantum codes",
    url="https://github.com/XanaduAI/ft-stack",
    packages=find_packages(),
    package_data={"ft_stack":["data/*"]},
    install_requires=[
        "certifi==2020.12.5",
        "cycler==0.10.0",
        "decorator==4.4.2",
        "kiwisolver==1.3.1",
        "matplotlib==3.3.3",
        "networkx==2.5",
        "numpy==1.21.0",
        "pandas==1.2.1",
        "Pillow==8.1.0",
        "pyparsing==2.4.7",
        "python-dateutil==2.8.1",
        "pytz==2020.5",
        "scipy==1.6.0",
        "six==1.15.0",
        "thewalrus==0.15.0",
        "pytest>=6.2.4",
        "pytest-cov>=2.12.1",
    ]
)
