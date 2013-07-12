#include "simple_stop_forwarder.hpp"

#include "timer.hpp"
#include "matrix.hpp"
#include "seq_io.hpp"
#include "debug.hpp"

#include <climits>
#include <vector>

namespace zipHMM {

  void SimpleStopForwarder::read_seq(const std::string &seq_filename, const size_t alphabet_size, std::vector<size_t> &nStatesSave) {
    orig_alphabet_size = alphabet_size;
    std::vector<unsigned> *prev_seq_p = 0;
    std::vector<unsigned> *current_seq_p = 0;
    std::vector<unsigned> *next_seq_p = 0;
    std::vector<size_t> pair_n;
    bool counted = false;
    size_t no_states_save = 0;
    size_t prev_max_count = UINT_MAX;
    bool nStatesSave_init_empty = nStatesSave.empty();
    if(!nStatesSave_init_empty)
      no_states_save = nStatesSave.back();
    
    prev_seq_p = new std::vector<unsigned>();
    zipHMM::readSeq(*prev_seq_p, seq_filename);
    orig_seq_length = prev_seq_p->size();

    current_seq_p = new std::vector<unsigned>();
    current_seq_p->assign(prev_seq_p->begin(), prev_seq_p->end()); // copy seq

    // weird special cases I need to handle
    if(orig_alphabet_size == 1 || no_states_save == 1 || current_seq_p->size() <= 2) {
      if(nStatesSave_init_empty) {
	nStates2seq[2] = *current_seq_p; // copy seq
	nStates2alphabet_size[2] = orig_alphabet_size;
	return;
      }
      for(std::vector<size_t>::iterator it = nStatesSave.begin(); it != nStatesSave.end(); ++it) {
	nStates2seq[*it] = *current_seq_p; // copy seq
	nStates2alphabet_size[*it] = orig_alphabet_size;
      }
      return;
    }

    // count each pair of symbols across the original sequence
    pair_n.resize(orig_alphabet_size * orig_alphabet_size);
    // first pair is special because we don't want to count overlapping pairs twice.
    update_pair_n(pair_n, s_pair( (*current_seq_p)[1] , (*current_seq_p)[2] ), orig_alphabet_size);
    counted = true;
    for(size_t i = 2; i < current_seq_p->size() - 1; ++i)
      add_count(pair_n, current_seq_p, i, counted, orig_alphabet_size);
    
    size_t new_alphabet_size = orig_alphabet_size;
    while(true) {
      // find pair with maximal number of non-overlapping occurrences
      s_pair max_pair;
      size_t max_count = 0;
      for(size_t i = 0; i < new_alphabet_size*new_alphabet_size; ++i) {
	if(pair_n[i] > max_count) {
	  max_count = pair_n[i];
	  unsigned left = unsigned(i / new_alphabet_size);
	  unsigned right = unsigned(i - left * new_alphabet_size);
	  max_pair = s_pair(left, right);
	}
      }

      // stopping criteria
      if(prev_max_count == max_count && nStatesSave_init_empty) {
	nStates2seq[2] = *current_seq_p; // copying sequence
	nStates2alphabet_size[2] = new_alphabet_size;
	break;
      } else if(max_count < no_states_save && !nStatesSave_init_empty) {
	nStates2seq[no_states_save] = *current_seq_p; // copying sequence
	nStates2alphabet_size[no_states_save] = new_alphabet_size;

	nStatesSave.pop_back();
	if(nStatesSave.empty())
	  break;
	else
	  no_states_save = nStatesSave.back();
      }

      // save the components of max_pair in symbol2pair
      symbol2pair[unsigned(new_alphabet_size)] = max_pair;
      
      pair_n.clear();
      pair_n.resize( (new_alphabet_size+1) * (new_alphabet_size+1));
      next_seq_p = new std::vector<unsigned>(current_seq_p->size() - max_count);
      (*next_seq_p)[0] = (*current_seq_p)[0]; // we never change the very first character
      // replace every occurrence of max_pair with a new alphabet symbol and count pairs at the same time
      size_t i = 1; // first position is special because we have not seen any pairs yet.
      size_t next_seq_pos = 1;
      if(matches(max_pair, current_seq_p, i)) { 
	(*next_seq_p)[next_seq_pos] = unsigned(new_alphabet_size);
	++i;
      } else {
	(*next_seq_p)[next_seq_pos] = (*current_seq_p)[i];
      }
      ++next_seq_pos;
      // do the rest of the sequence
      for(i = i+1; i < current_seq_p->size() - 1; i++) {
	if(matches(max_pair, current_seq_p, i)) {
	(*next_seq_p)[next_seq_pos] = unsigned(new_alphabet_size);
	  ++i;
	} else {
	(*next_seq_p)[next_seq_pos] = (*current_seq_p)[i];
	}

	add_count(pair_n, next_seq_p, next_seq_pos - 1, counted, new_alphabet_size + 1);
	++next_seq_pos;
      }
      if(i == current_seq_p->size() - 1) { // add last symbol of sequence
	(*next_seq_p)[next_seq_pos] = (*current_seq_p)[i];
	add_count(pair_n, next_seq_p, next_seq_pos - 1, counted, new_alphabet_size + 1);
      }
	    
      // set up for next round
      prev_max_count = max_count;
      delete prev_seq_p;
      prev_seq_p = current_seq_p;
      current_seq_p = next_seq_p;
      new_alphabet_size++;
    }

    delete prev_seq_p;
    delete current_seq_p;
    // delete next_seq_p; has already been freed!
  }

  void SimpleStopForwarder::read_seq(const std::string &seq_filename, const size_t alphabet_size, const size_t no_states) {
    std::vector<size_t> nStatesSave;
    nStatesSave.push_back(no_states);
    read_seq(seq_filename, alphabet_size, nStatesSave);
  }

  void SimpleStopForwarder::read_seq(const std::string &seq_filename, const size_t alphabet_size) {
    std::vector<size_t> nStatesSave;
    read_seq(seq_filename, alphabet_size, nStatesSave);
  }




} // namespace
