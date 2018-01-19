import os
import glob
import numpy
from setuptools import setup, Extension
import pysam

src_list = glob.glob("src/*.c")
pyx_list = glob.glob("dnaUtilsPy/*.pyx")
include_dirs = pysam.get_include() + ["src", numpy.get_include()]
extra_link_args = pysam.get_libraries()

extensions = []
for pyx in pyx_list:
    assert pyx.endswith(".pyx")
    submodule = os.path.basename(pyx)[:-len(".pyx")]
    extensions.append(
        Extension(
            "dnaUtilsPy." + submodule,
            sources=[pyx] + src_list,
            include_dirs=include_dirs,
            extra_link_args=extra_link_args))

setup(
    name="dnaUtilsPy",
    packages=["dnaUtilsPy"],
    ext_modules=extensions,
    install_requires=['numpy', 'pysam'],
    setup_requires=['numpy', 'pysam', 'cython'])
