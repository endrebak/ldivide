import os
import sys

from distutils.core import setup

from setuptools import find_packages, Extension, Command
from Cython.Build import cythonize

__version__ = open("ldivide/version.py").readline().split(" = ")[1].replace(
    '"', '').strip()
macros = []

setup_requires = ["cython"]
install_requires = [
    "numpy"
]


compile_options = [
    "-Ofast", "-Wall"
]  #, "-frename-registers", "-funroll-loops"] # , "-lgzstream", "-lz"

# print(conda_lib)

extensions = [
    Extension(
        "ldivide.src.calc_covar", ["ldivide/src/calc_covar.pyx"],
        language="c"),
    Extension(
        "ldivide.src.calc_autocovar", ["ldivide/src/calc_autocovar.pyx"],
        language="c"),
    Extension(
        "ldivide.src.calc_square", ["ldivide/src/calc_square.pyx"],
        language="c"),
    # Extension(
    #     "ldetect2.src.matrix_to_vector", ["ldetect2/src/matrix_to_vector.pyx"],
    #     language="c"),
    # Extension(
    #     "ldetect2.src.m2v", ["ldetect2/src/m2v.pyx"],
    #     language="c"),
]


setup(
    name="ldivide",
    packages=find_packages(),
    ext_modules=cythonize(extensions, annotate=True, language_level='3'),
    version=__version__,
    description="Fast finding of approx independent linkage disequilibrium blocks in human populations.",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/ldivide",
    keywords=["genetics"],
    license=["BSD"],
    setup_requires=setup_requires,
    install_requires=install_requires,
    include_dirs=["."],
    long_description="."
    )
