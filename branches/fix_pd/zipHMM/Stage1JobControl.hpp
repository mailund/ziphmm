#ifndef STAGE_1_JOB_CONTROL_HPP
#define STAGE_1_JOB_CONTROL_HPP

#include "debug.hpp"
#include "hmm_utils.hpp"

#include <map>
#include <deque>
#include <pthread.h>
#include <stdint.h>
#include <stack>

namespace zipHMM {

  struct Stage1JobControl {
    bool *symbol2is_computed;
    bool *symbol2is_added;
    const std::map<unsigned, s_pair> *symbol2pair;
    std::map<unsigned, std::deque<unsigned> > left2pair_symbols;
    std::map<unsigned, std::deque<unsigned> > right2pair_symbols;
    std::stack<unsigned> compute_stack;
    size_t orig_alphabet_size;
    size_t alphabet_size;
    size_t no_fetched_symbols;
    size_t no_threads;

    pthread_mutex_t fetch_lock;
    pthread_mutex_t update_lock;
    pthread_cond_t empty_stack_cond;

    Stage1JobControl(const size_t orig_alphabet_size, const size_t alphabet_size, const std::map<unsigned, s_pair> *symbol2pair, size_t no_threads);

    ~Stage1JobControl() {
      pthread_mutex_destroy(&fetch_lock);
      pthread_mutex_destroy(&update_lock);
      pthread_cond_destroy(&empty_stack_cond);

      delete[] symbol2is_computed;
      delete[] symbol2is_added;
    }
    
    void inline lock_fetch()    { pthread_mutex_lock  (&fetch_lock);  }
    void inline unlock_fetch()  { pthread_mutex_unlock(&fetch_lock);  }
    void inline lock_update()   { pthread_mutex_lock  (&update_lock); }
    void inline unlock_update() { pthread_mutex_unlock(&update_lock); }

    void inline wait_empty_stack()   { pthread_cond_wait(&empty_stack_cond, &fetch_lock); }
    void inline signal_empty_stack() { pthread_cond_signal(&empty_stack_cond); }
    
    bool inline &is_computed(unsigned symbol) const { return symbol2is_computed[symbol]; }
    void inline set_computed(unsigned symbol, bool val) { symbol2is_computed[symbol] = val; }
    const size_t inline getOriginalAlphabetSize() const { return orig_alphabet_size; }
    const size_t inline getAlphabetSize() const { return alphabet_size; }
    const std::deque<unsigned> inline &get_pairs_from_left(unsigned symbol) const {  return left2pair_symbols.find(symbol)->second; }
    const std::deque<unsigned> inline &get_pairs_from_right(unsigned symbol) const { return right2pair_symbols.find(symbol)->second; }
    
    bool inline is_done() const { return (no_fetched_symbols == alphabet_size); }

    size_t get_jobs(std::deque<unsigned> &result);

    void update_symbol2is_computed(const std::deque<unsigned> &jobs_done);

    void find_jobs_to_add(const std::deque<unsigned> &jobs_done, std::deque<unsigned> &jobs_to_add);

    void add_jobs(const std::deque<unsigned> &jobs_to_add);
    
    void update(const std::deque<unsigned> &jobs_done) {
      update_symbol2is_computed(jobs_done);

      std::deque<unsigned> jobs_to_add;
      find_jobs_to_add(jobs_done, jobs_to_add);

      add_jobs(jobs_to_add);
    }

  }; // struct

}

#endif
