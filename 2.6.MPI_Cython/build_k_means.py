from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext

# https://software.intel.com/en-us/articles/thread-parallelism-in-cython
setup(
    name = "k_means_cython",
    cmdclass = {"build_ext": build_ext},
    ext_modules =
        [
            Extension("k_means_cython",
                ["k_means_cython.pyx"],
                extra_compile_args = ["-O3", "-ipo", "-fopenmp"],
                extra_link_args=['-fopenmp']
        )
    ]
)

