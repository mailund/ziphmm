#include "prob_spaces.hpp"

#include <limits>

namespace zipHMM {
  const double LinearSpace::ZERO = 0;
  const double LinearSpace::ONE = 1;
  
  const double LogSpace::ZERO = -std::numeric_limits<double>::infinity();
  const double LogSpace::ONE = 0;
}
