#include "hmm_utils.hpp"
#include "matrix.hpp"

#include <vector>

namespace zipHMM {

  bool matches(const s_pair &p, const std::vector<unsigned> *seq_p, const size_t idx) {
    return p.first == (*seq_p)[idx] && p.second == (*seq_p)[idx + 1];
  }

  void update_pair_n(std::vector<size_t> &pair_n, const s_pair &pair, const size_t alphabet_size) {
    unsigned left  = pair.first;
    unsigned right = pair.second;
    size_t idx = left * alphabet_size + right;
    pair_n[idx]++;
  }

  /**
   *  Two overlapping occurrences of the same pair of characters
   *  cannot both be replaced by the same new character. We therefore
   *  add one to pair_n[ (seq[position], seq[position+1]) ] iff an
   *  identical pair was not counted in the previous position. That is
   *  either if this is the first occurrence of this pair in the
   *  sequence, OR if it does not match in the previous position, OR
   *  if it was not counted in the previous position (due to this
   *  restriction).
   */
  void add_count(std::vector<size_t> &pair_n, std::vector<unsigned> *seq_p, size_t position, bool &counted, size_t alphabet_size) {
    s_pair current_pair((*seq_p)[position], (*seq_p)[position+1]);
    // if it didn't match previous position, increase count
    if(!matches(current_pair, seq_p, position-1)) { 
      update_pair_n(pair_n, current_pair, alphabet_size);
      counted = true;
      // if it matched at the previous position, check if the previous position was counted before increasing count
    } else if(!counted) { 
      update_pair_n(pair_n, current_pair, alphabet_size);
      counted = true;
      // if it did match and was counted at previous position, don't increase count at this position.
    } else {
      counted = false;
    }
  }

  double init_apply_em_prob(Matrix &res, const Matrix &pi, const Matrix &B, const unsigned symbol) {
    const size_t nStates = pi.get_height();
    res.reset(nStates, 1);
    
    for(size_t i = 0; i < nStates; ++i)
      res(i, 0) = pi(i, 0) * B(i, symbol);
    
    return res.normalize();
  }
  
  double apply_em_prob(Matrix &res, const Matrix &A, const Matrix &B, const unsigned symbol) {
    const size_t nStates = A.get_height();
    res.reset(nStates, nStates);
    
    for(size_t r = 0; r < nStates; ++r) {
      for(size_t c = 0; c < nStates; ++c)
	res(r, c) = A(c, r) * B(r, symbol);
    }

    return res.normalize();
  }

  void make_em_trans_probs_array(double *symbol2scale, Matrix *symbol2matrix, const Matrix &A, const Matrix &B) {
    for(size_t i = 0; i < B.get_width(); ++i) {
      symbol2scale[unsigned(i)] = std::log( apply_em_prob(symbol2matrix[unsigned(i)], A, B, unsigned(i)) );
    }
  }

}
