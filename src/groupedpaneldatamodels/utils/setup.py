from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

ext_modules = [
    Extension(
        name="bonhomme_manresa",
        sources=["bonhomme_manresa.pyx"],
        include_dirs=[np.get_include()],
        define_macros=[("CYTHON_TRACE", "1")],
        extra_compile_args=["-O2"],
    )
]

setup(
    name="bonhomme_manresa",
    ext_modules=cythonize(
        ext_modules,
        language_level=3,
        compiler_directives={
            "boundscheck": False,
            "wraparound": False,
            "nonecheck": False,
            "cdivision": True,
        },
    ),
    zip_safe=False,
)
