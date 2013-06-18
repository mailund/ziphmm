#include "test_utils.hpp"

#include "performance_description.hpp"

#include <sstream>
#include <fstream>
#include <iostream>

namespace zipHMM {
    const double EPS = 0.0000000001;

  DeviceDescriptor createDeviceDescriptor(size_t nDevices) {
    return DeviceDescriptor(PTHREAD_PROCESSING_DEVICE_FACTORY, nDevices);
  }

  void assertTrue(bool val, const std::string &text) {
    if(!val) {
      std::cerr << text << std::endl;
      std::exit(-1);
    }
  }
  
  void assertClose(double x1, double x2, const std::string &text, double eps) {
    double maxDif = std::max(eps, std::min(fabs(x1), fabs(x2)) * eps);
    if(fabs(x2 - x1) > maxDif) {
      std::cerr << text << " Found " << x1 << " but expected " << x2 << std::endl;
      std::exit(-1);
    }
  }

  void assertClose(const Matrix &found, const Matrix &expected, double eps) {
    assertTrue(found.get_height() == expected.get_height(), "Matrix height differ.");
    assertTrue(found.get_width() == expected.get_width(), "Matrix width differ.");
    for(unsigned i = 0; i < found.get_height(); ++i) {
      for(unsigned j = 0; j < found.get_width(); ++j)
	assertClose(found(i, j), expected(i, j), "Matrix entries differ.", eps);
    }
  }
  
  void assertEqual(const std::vector<unsigned> &vec, const std::string &str) {
    std::ostringstream vecStr;
    std::string sep = "";
    for(unsigned i = 0; i < vec.size(); ++i) {
      vecStr << sep << vec[i];
      sep = " ";
    }
    
    assertEqual<std::string>(vecStr.str(), str);
  }

    void assertFilesEqual(const std::string &fn1, const std::string &fn2) {
      std::ifstream i1(fn1.c_str()), i2(fn2.c_str());
      
      while(i1 && i2) {
	if(i1.get() != i2.get()) {
	  std::cerr << "Read and written files are not equal!" << std::endl;
	  std::exit(-1);
	}
      }

      if(i1.good() != i2.good()) {
	std::cerr << "Read and written files did not end at the same time!" << std::endl;
	std::exit(-1);
      }
      
      i1.close();
      i2.close();
    }
}
