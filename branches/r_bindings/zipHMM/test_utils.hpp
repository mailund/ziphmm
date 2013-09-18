#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include "matrix.hpp"
#include "performance_description.hpp"

#include <cstdlib>

namespace zipHMM {
  extern const double EPS;
  
  DeviceDescriptor createDeviceDescriptor(unsigned nDevices);
  
  void assertTrue(bool val, const std::string &text = "Assertion failed!");
  void assertClose(double x1, double x2, const std::string &text = "Assertion failed:", double eps = EPS);
  void assertClose(const Matrix &found, const Matrix &expected, double eps = EPS);
  void assertEqual(const std::vector<unsigned> &vec, const std::string &str);
  void assertFilesEqual(const std::string &fn1, const std::string &fn2);
  
  template<typename T>
  void assertEqual(const T& x1, const T& x2, const std::string &text = "Assertion failed:") {
    if(x1 != x2) {
      std::cerr << text << " Found " << x1 << " but expected " << x2 << std::endl;
      std::exit(-1);
    }
  }
}

#endif
