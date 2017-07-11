#!/usr/bin/env python

from distutils.core import setup,Command
from numpy.distutils.misc_util import get_numpy_include_dirs
from distutils.extension import Extension
import os

# you can customize
#os.environ["CC"] = "clang"
#os.environ["CXX"] = "clang"

# cython extension for the euler integration
eulerint = Extension("eulerint",
                    include_dirs = get_numpy_include_dirs(),
                    sources = ["eulerint.c"])

class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        import sys,subprocess
        errno = subprocess.call([sys.executable,'tests/runtests.py'])
        raise SystemExit(errno)

setup(name='pydynet',
      version='1.3.1',
      description='Python Package for Dynamics on Networks',
      author='Kevin Brown and Ann Hermundstad',
      author_email='kevin.s.brown@uconn.edu',
      url='https://github.com/thelahunginjeet/pydynet',
      packages=['pydynet'],
      package_dir = {'pydynet': ''},
      package_data = {'pydynet' : ['tests/*.py']},
      ext_package = 'pydynet',
      ext_modules = [eulerint],
      cmdclass = {'test': PyTest},
      license='BSD-3',
      classifiers=[
          'License :: OSI Approved :: BSD-3 License',
          'Intended Audience :: Developers',
          'Intended Audience :: Scientists',
          'Programming Language :: Python',
          'Topic :: Graph Theory',
          'Topic :: Dynamical Systems'
      ],
    )
