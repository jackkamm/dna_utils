import os
import numpy
from setuptools import setup, Extension
import pysam

extensions = [
    Extension(
        "dnaUtilsPy.vcf_pairwise_diffs.wrapper",
        sources=["dnaUtilsPy/vcf_pairwise_diffs/wrapper.pyx",
                 "src/vcfpairdiffs.c"],
        include_dirs=pysam.get_include() + [
            "src", numpy.get_include()],
        extra_link_args=pysam.get_libraries()),
]

setup(name="dnaUtilsPy",
	  packages=["dnaUtilsPy", "dnaUtilsPy.vcf_pairwise_diffs"],
      ext_modules=extensions,
      install_requires=['numpy', 'pysam'],
      setup_requires=['numpy', 'pysam', 'cython'])
