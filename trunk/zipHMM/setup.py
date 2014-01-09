from distutils.core import setup
from platform import system

setup(name = "zipHMM",
      version = "1.0.2",
      description = "zipHMM is a package for working with general discrete hidden Markov models, greatly speeding up the likelihood computation by using the zipForward algorithm.",
      author = "Andreas Sand",
      author_email = "asand@birc.au.dk",
      url = "http://birc.au.dk/software/zipHMM",
      packages = ["pyZipHMM"],
      package_data = {"": ["libpyZipHMM.so"]},
      license = "LGPL" 
    )
