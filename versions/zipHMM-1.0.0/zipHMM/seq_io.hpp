#ifndef SEQ_IO_HPP
#define SEQ_IO_HPP

#include <string>
#include <vector>
#include <iostream>

namespace zipHMM {
  void readSeq(std::vector<unsigned> &result, const std::string &filename);

  void readSeq(std::vector<unsigned> &result, std::istream &in);

  void writeSeq(const std::vector<unsigned> &seq, const std::string &filename);

  void writeSeq(const std::vector<unsigned> &seq, std::ostream &out);
}

#endif
