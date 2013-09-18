#ifndef SIMPLE_STOP_FORWARDER_HPP
#define SIMPLE_STOP_FORWARDER_HPP

#include "forwarder.hpp"

#include <vector>
#include <string>

namespace zipHMM {

  class SimpleStopForwarder : public Forwarder {

  public:

    SimpleStopForwarder() : Forwarder() { }

    void read_seq(const std::string &seq_filename, const size_t alphabet_size);
    void read_seq(const std::string &seq_filename, const size_t alphabet_size, std::vector<size_t> &nStatesSave);
    void read_seq(const std::string &seq_filename, const size_t alphabet_size, const size_t no_states);

    ~SimpleStopForwarder() { }

  };

} // namespace

#endif
