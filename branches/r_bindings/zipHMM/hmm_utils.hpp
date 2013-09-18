#ifndef HMM_UTILS_HPP
#define HMM_UTILS_HPP

#include "matrix.hpp"

#include <map>
#include <vector>

namespace zipHMM {

  typedef std::pair<unsigned, unsigned> s_pair;
  
  bool matches(const s_pair &pair, const std::vector<unsigned> *seq, const size_t idx);

  void add_count(std::vector<size_t> &pair_n, std::vector<unsigned> *seq, size_t position, bool &counted, size_t alphabet_size);
  void update_pair_n(std::vector<size_t> &pair_n, const s_pair &pair, const size_t alphabet_size);
  
  double init_apply_em_prob(Matrix &res, const Matrix &pi, const Matrix &B, const unsigned symbol);
  double apply_em_prob(Matrix &res, const Matrix &A, const Matrix &B, const unsigned symbol);
  void make_em_trans_probs_array(double *symbol2scale, Matrix *symbol2matrix, const Matrix &A, const Matrix &B);
}

#endif
