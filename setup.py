#!/usr/bin/env python
# coding: utf-8
# MIT License

"""
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.

---
This code borrows heavily from 
https://github.com/qmlcode/qml

License: 
---

#
# Copyright (c) 2018 Lars Andersen Bratholm
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""

import sys
from numpy.distutils.core import Extension, setup

from mkldiscover import mkl_exists

FORTRAN = "f90"

# GNU (default)
COMPILER_FLAGS = ["-O3", "-fopenmp", "-m64", "-march=native", "-fPIC",
                    "-Wno-maybe-uninitialized", "-Wno-unused-function", "-Wno-cpp"]
LINKER_FLAGS = ["-lgomp"]
# MATH_LINKER_FLAGS = ["-lblas", "-llapack", "-latlas", "-fopenmp"]
MATH_LINKER_FLAGS = ["-lblas", "-llapack"]

# UNCOMMENT TO FORCE LINKING TO MKL with GNU compilers:
if mkl_exists(verbose=True):
    LINKER_FLAGS = ["-lgomp", "-lpthread", "-lm", "-ldl"]
    MATH_LINKER_FLAGS = ["-L${MKLROOT}/lib/intel64", "-lmkl_rt"]

# For clang without OpenMP: (i.e. most Apple/mac system)
if sys.platform == "darwin" and all(["gnu" not in arg for arg in sys.argv]):
    COMPILER_FLAGS = ["-O3", "-m64", "-march=native", "-fPIC"]
    LINKER_FLAGS = []
    # MATH_LINKER_FLAGS = ["-lblas", "-llapack", "-latlas", "-fopenmp"]
    MATH_LINKER_FLAGS = ["-lblas", "-llapack"]


# Intel
if any(["intelem" in arg for arg in sys.argv]):
    COMPILER_FLAGS = ["-xHost", "-O3", "-axAVX", "-qopenmp"]
    LINKER_FLAGS = ["-liomp5", "-lpthread", "-lm", "-ldl"]
    MATH_LINKER_FLAGS = ["-L${MKLROOT}/lib/intel64", "-lmkl_rt"]


ext_fgpfchl = Extension(name = '.fgpfchl.frepresentations',
                          sources = [
                          'gplearn/fgpfchl/frepresentations.f90',
                              ],
                          extra_f90_compile_args = COMPILER_FLAGS,
                          extra_f77_compile_args = COMPILER_FLAGS,
                          extra_compile_args = COMPILER_FLAGS,
                          extra_link_args = LINKER_FLAGS + MATH_LINKER_FLAGS,
                          language = FORTRAN,
                          f2py_options=['--quiet'])

def requirements():
    with open('requirements.txt') as f:
        return [line.rstrip() for line in f]


# use README.md as long description
def readme():
    with open('README.rst') as f:
        return f.read()


def setup_gplearn():

    setup(

        name="gplearn",
        packages=[
            'gplearn',
            ],

        # metadata
        version="0.0.1",
        platforms = 'Any',
        keywords = ['Machine Learning', 'Group Period'],
        classifiers = [],

        # set up package contents

        ext_package = 'gplearn',
        ext_modules = [
              ext_fgpfchl
        ],
)


if __name__ == '__main__':

    setup_gplearn()