#include "Stage1JobControl.hpp"
#include "debug.hpp"
#include "hmm_utils.hpp"

#include <map>
#include <deque>
#include <pthread.h>
#include <stdint.h>
#include <stack>

namespace zipHMM {

  Stage1JobControl::Stage1JobControl(const size_t orig_alphabet_size, const size_t alphabet_size, const std::map<unsigned, s_pair> *symbol2pair, size_t no_threads)
    : symbol2pair(symbol2pair),
      orig_alphabet_size(orig_alphabet_size), alphabet_size(alphabet_size), 
      no_fetched_symbols(orig_alphabet_size), no_threads(no_threads) {
      
    symbol2is_computed = new bool[alphabet_size];
    symbol2is_added = new bool[alphabet_size];

    for(size_t i = 0; i < orig_alphabet_size; ++i) {
      set_computed((unsigned) i, true);
      symbol2is_added[i] = true;
      left2pair_symbols[(unsigned) i]; //initialize key with empty deque
      right2pair_symbols[(unsigned) i]; //initialize key with empty deque
    }

    for(size_t i = orig_alphabet_size; i < alphabet_size; ++i) {
      set_computed((unsigned) i, false);
      symbol2is_added[i] = false;
      left2pair_symbols[(unsigned) i]; // initialize key with empty deque
      right2pair_symbols[(unsigned) i];  // initialize key with empty deque
    }

    // initialize compute stack data structures
    for(size_t i = orig_alphabet_size; i < alphabet_size; ++i) {
      s_pair symbol_pair = symbol2pair->find((unsigned) i)->second;
      unsigned left  = symbol_pair.second;
      unsigned right = symbol_pair.first;
	
      if(symbol2is_computed[left] && symbol2is_computed[right]) {
	compute_stack.push((unsigned) i);	
	symbol2is_added[i] = true;
      }
	
      left2pair_symbols[left].push_back((unsigned) i);
      right2pair_symbols[right].push_back((unsigned) i);
    }

    pthread_mutex_init(&fetch_lock, 0);
    pthread_mutex_init(&update_lock, 0);
    pthread_cond_init(&empty_stack_cond, 0);
  }

  size_t Stage1JobControl::get_jobs(std::deque<unsigned> &result) {
    lock_fetch();
    while(compute_stack.empty() && !is_done()) {
      wait_empty_stack();
    }
    if(is_done()) {
      signal_empty_stack();
      unlock_fetch();
      return 0;
    }

    size_t no_symbols_to_fetch = compute_stack.size() / no_threads;
    if(compute_stack.size() % no_threads != 0)
      no_symbols_to_fetch++;
      
    lock_update();
    for(size_t i = 0; i < no_symbols_to_fetch; ++i) {
      result.push_back(compute_stack.top());
      compute_stack.pop();
    }
    unlock_update();

    no_fetched_symbols += no_symbols_to_fetch;
    DEBUG("alphabet_size: %d\n", alphabet_size);
    DEBUG("fetched %d symbols. no_fetched_symbols now %d : %d\n", no_symbols_to_fetch, no_fetched_symbols, alphabet_size);

    unlock_fetch();
    return no_symbols_to_fetch;
  }

  void Stage1JobControl::update_symbol2is_computed(const std::deque<unsigned> &jobs_done) {
    lock_update();
    for(std::deque<unsigned>::const_iterator it = jobs_done.begin(); it != jobs_done.end(); ++it)
      set_computed(*it, true);
    unlock_update();
    DEBUG("-- : %d symbols completed\n", jobs_done.size());
  }

  void Stage1JobControl::find_jobs_to_add(const std::deque<unsigned> &jobs_done, std::deque<unsigned> &jobs_to_add) {
    for(std::deque<unsigned>::const_iterator it = jobs_done.begin(); it != jobs_done.end(); ++it) {
      unsigned symbol = (*it);
	
      for(std::deque<unsigned>::const_iterator it = get_pairs_from_left(symbol).begin(); it != get_pairs_from_left(symbol).end(); ++it) {
	s_pair left_symbol_pair = symbol2pair->find(*it)->second;
	if(std::find(jobs_to_add.begin(), jobs_to_add.end(), *it) == jobs_to_add.end() && is_computed(left_symbol_pair.first) ) { 
	  jobs_to_add.push_back(*it);
	  DEBUG("left suggested %d\n", *it);
	}
      }
	
      for(std::deque<unsigned>::const_iterator it = get_pairs_from_right(symbol).begin(); it != get_pairs_from_right(symbol).end(); ++it) {
	s_pair right_symbol_pair = symbol2pair->find(*it)->second;
	if(std::find(jobs_to_add.begin(), jobs_to_add.end(), *it) == jobs_to_add.end() && is_computed(right_symbol_pair.second) ) {
	  jobs_to_add.push_back(*it);
	  DEBUG("right suggested %d\n", *it);
	}
      }
    }
  }

  void Stage1JobControl::add_jobs(const std::deque<unsigned> &jobs_to_add) {
    lock_update();
    for(std::deque<unsigned>::const_iterator it = jobs_to_add.begin(); it != jobs_to_add.end(); ++it) {
      if(!symbol2is_added[*it]) {
	compute_stack.push(*it);
	symbol2is_added[*it] = true;
	DEBUG("pushed %d\n", *it);
      }
    }
    DEBUG("%d symbols pushed\n", jobs_to_add.size());
    signal_empty_stack();
    unlock_update();
  }
    
}
