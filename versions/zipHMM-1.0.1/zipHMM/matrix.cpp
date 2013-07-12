#include "matrix.hpp"

std::ostream &operator<<(std::ostream &out, const zipHMM::Matrix &mat) {
  for(unsigned i = 0; i < mat.get_height(); ++i) {
    for(unsigned j = 0; j < mat.get_width(); ++j)
      out << mat(i, j) << " ";
    out << std::endl;
  }
  
  return out;
}

