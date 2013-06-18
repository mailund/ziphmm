#include "seq_io.hpp"

#include <fstream>
#include <cstdlib>
#include <vector>

namespace zipHMM {

    void readSeq(std::vector<unsigned> &result, const std::string &filename) {
      std::ifstream in(filename.c_str());
      
      if(!in) {
	std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
	exit(-1);
      }
      
      readSeq(result, in);
      
      in.close();
    }

    void readSeq(std::vector<unsigned> &result, std::istream &in) {
      result.clear();
      
      unsigned tmp = 0;
      while(in >> tmp)
	result.push_back(tmp);
    }

    void writeSeq(const std::vector<unsigned> &seq, const std::string &filename) {
      std::ofstream out(filename.c_str());
      
      if(!out) {
	std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
	exit(-1);
      }

      writeSeq(seq, out);
      
      out.close();
    }

    void writeSeq(const std::vector<unsigned> &seq, std::ostream &out) {
      std::string sep = "";
      for(unsigned i = 0; i < seq.size(); ++i) {
	out << sep << seq[i];
	sep = " ";
      }
      out << std::endl;
    }
}
