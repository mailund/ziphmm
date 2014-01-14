from distutils.core import setup
from distutils.command.build import build
from platform import system
from os.path import expanduser

class calibrateRunner(build):

    def run(self):
        """Calibrating pyZipHMM for multi-threaded applications"""
        from pyZipHMM.pyZipHMM import calibrate

        print "Calibrating pyZipHMM for multi-threaded applications."
        filePath = expanduser("~") + "/.ziphmm.devices"
        calibrate(filePath)
        print "Calibration saved to preference file: %s" % filePath

        build.run(self)
        

setup(name = "zipHMM",
      version = "1.0.2",
      description = "zipHMM is a package for working with general discrete hidden Markov models, greatly speeding up the likelihood computation by using the zipForward algorithm.",
      author = "Andreas Sand",
      author_email = "asand@birc.au.dk",
      url = "http://birc.au.dk/software/zipHMM",
      packages = ["pyZipHMM"],
      package_data = { "": ["libpyZipHMM.so"]},
      license = "LGPL",
      cmdclass = { 'build' : calibrateRunner }
    )
