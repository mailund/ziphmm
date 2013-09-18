#ifndef PROB_SPACES_HPP
#define PROB_SPACES_HPP

#include <cmath>
#include <algorithm>

namespace zipHMM {

    class LinearSpace {
    public:
      static const double ZERO;
      static const double ONE;
      static double fromLinear(double x) { return x; };
      static double toLogSpace(double x) { return log(x); };
      static double add(double a, double b) { return a + b; };
      static double mult(double a, double b) { return a * b; };
      static double div(double a, double b) { return a / b; };
      static double max(double a, double b) { return std::max(a, b); };
    };
  
  class LogSpace {
  public:
    static const double ZERO;
    static const double ONE;
    static double fromLinear(double x) { return log(x); };
    static double toLogSpace(double x) { return x; };
    static double add(double a, double b) {
      if(a == ZERO)
	return b;
      else if(b == ZERO)
	return a;
      else if(a > b)
	return a + log1p(exp(b - a));
      else
	return b + log1p(exp(a - b));
    };
    static double mult(double a, double b) { return a + b; };
    static double div(double a, double b) { return a - b; };
    static double max(double a, double b) { return std::max(a, b); };
  };
}

#endif
