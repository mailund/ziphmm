#define COMPLEX COMPLEX1 // COMPLEX is declared in two third-party header files. This gets rid of one of them

#include "forwarder.hpp"
#include "simple_stop_forwarder.hpp"
#include "simple_forwarder.hpp"
#include "matrix.hpp"
#include "hmm_io.hpp"
#include "posterior_decoding.hpp"
#include "viterbi.hpp"
#include "calibrate.hpp"

#undef COMPLEX

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <cstdlib>
#include <vector>

size_t &parse_integer(SEXP sexp, std::string name) {
  if(!isInteger(sexp)) {
    std::string error_msg = name + " should be an integer.";
    error(error_msg.c_str());
  }
  size_t *int_ptr = (size_t *) INTEGER(sexp);
  return *int_ptr;
}

std::string parse_string(SEXP sexp, std::string name) {
  if(!isString(sexp) || length(sexp) < 1) {
    std::string error_msg = name + "should be a string.";
    error(error_msg.c_str());
  }
  std::string str = CHAR(STRING_ELT(sexp,0));
  return str;
}

void parse_HMM_parameters(SEXP pi_sexp, SEXP A_sexp, SEXP B_sexp, double **pi_content, double **A_content, double **B_content) {
  if(!isVector(pi_sexp))
    error("pi should be a vector");
  *pi_content = REAL(pi_sexp);
  if(!isMatrix(A_sexp))
    error("A should be a matrix");
  *A_content = REAL(A_sexp);
  if(!isMatrix(B_sexp))
    error("B should be a matrix");
  *B_content = REAL(B_sexp);
}

void convert_HMM_parameters(SEXP pi_sexp, SEXP A_sexp, SEXP B_sexp, SEXP no_states_sexp, SEXP alphabet_size_sexp,
	       zipHMM::Matrix &pi, zipHMM::Matrix &A, zipHMM::Matrix &B) {
  int no_states = parse_integer(no_states_sexp, "no_states");
  int alphabet_size = parse_integer(alphabet_size_sexp, "alphabet_size");
  
  double *pi_content, *A_content, *B_content;
  parse_HMM_parameters(pi_sexp, A_sexp, B_sexp, &pi_content, &A_content, &B_content);

  pi.reset(no_states, 1);
  A.reset(no_states, no_states);
  B.reset(no_states, alphabet_size);
  
  for(size_t r = 0; r < no_states; ++r) {
    pi(r, 0) = pi_content[r];
   
    for(size_t c = 0; c < no_states; ++c)
      A(r,c) = A_content[c*no_states + r];
    
    for(size_t c = 0; c < alphabet_size; ++c)
      B(r,c) = B_content[c*no_states + r];
  }
}

// HMM_IO
extern "C" {
  void c_read_HMM_spec(SEXP no_states_sexp, SEXP alphabet_size_sexp, SEXP filename_sexp) {

    size_t &no_states = parse_integer(no_states_sexp, "no_states");
    size_t &alphabet_size = parse_integer(alphabet_size_sexp, "alphabet_size");
    std::string filename = parse_string(filename_sexp, "filename");

    size_t c_no_states;
    size_t c_alphabet_size;
    zipHMM::read_HMM_spec(c_no_states, c_alphabet_size, std::string(filename));

    no_states = c_no_states;
    alphabet_size = c_alphabet_size;
  }

  void c_read_HMM(SEXP pi_sexp, SEXP A_sexp, SEXP B_sexp, SEXP filename_sexp) {
    double *pi_content, *A_content, *B_content;
    parse_HMM_parameters(pi_sexp, A_sexp, B_sexp, &pi_content, &A_content, &B_content);

    std::string filename = parse_string(filename_sexp, "filename");

    zipHMM::Matrix c_pi, c_A, c_B;
    zipHMM::read_HMM(c_pi, c_A, c_B, std::string(filename));
    
    for(size_t i = 0; i < c_pi.get_height(); ++i)
      pi_content[i] = c_pi(i,0);

    for(size_t r = 0; r < c_A.get_height(); ++r)
      for(size_t c = 0; c < c_A.get_width(); ++c)
    	A_content[c * c_A.get_height() + r] = c_A(r,c);

    for(size_t r = 0; r < c_B.get_height(); ++r)
      for(size_t c = 0; c < c_B.get_width(); ++c)
    	B_content[c * c_B.get_height() + r] = c_B(r,c);
  }

  void c_write_HMM(SEXP pi_sexp, SEXP A_sexp, SEXP B_sexp, SEXP no_states_sexp, SEXP alphabet_size_sexp, SEXP filename_sexp) {
    double *pi_content, *A_content, *B_content;
    parse_HMM_parameters(pi_sexp, A_sexp, B_sexp, &pi_content, &A_content, &B_content);

    std::string filename = parse_string(filename_sexp, "filename");

    zipHMM::Matrix c_pi, c_A, c_B;
    convert_HMM_parameters(pi_sexp, A_sexp, B_sexp, no_states_sexp, alphabet_size_sexp, c_pi, c_A, c_B);
    
    write_HMM(c_pi, c_A, c_B, filename);
  }

} // end HMM_IO

zipHMM::Forwarder *cast_forwarder(SEXP f_ptr_sexp) {
    void *f_cptr = R_ExternalPtrAddr(f_ptr_sexp);
    zipHMM::Forwarder *f = static_cast<zipHMM::Forwarder *>(f_cptr);
    
    return f;
}


// Forwarder
extern "C" {

  SEXP Forwarder_new() {
    zipHMM::Forwarder *f = new zipHMM::Forwarder();
    return R_MakeExternalPtr((void *)f, R_NilValue, R_NilValue);
  }

  void Forwarder_read_seq(SEXP f_ptr_sexp, SEXP seq_filename_sexp, SEXP alphabet_size_sexp, SEXP nStatesSave_sexp, SEXP min_no_evals_sexp) {
    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);

    std::string seq_filename = parse_string(seq_filename_sexp, "seq_filename");
    int alphabet_size = parse_integer(alphabet_size_sexp, "alphabet_size");
    int min_no_evals = parse_integer(min_no_evals_sexp, "min_no_evals");

    std::vector<size_t> nStatesSave;
    if(!isNull(nStatesSave_sexp)) {
      if(isVector(nStatesSave_sexp)) {
	int *nStatesSave_content = INTEGER(nStatesSave_sexp);
	for(size_t i = 0; i < GET_LENGTH(nStatesSave_sexp); ++i) // something like this
	  nStatesSave.push_back(nStatesSave_content[i]);
      } else {
	error("nStatesSave_sexp should be a vector");
      }
    }
    
    f->read_seq(seq_filename, alphabet_size, nStatesSave, min_no_evals);
  }

  void Forwarder_read_from_directory(SEXP f_ptr_sexp, SEXP directory_sexp) {
    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);
    std::string directory = parse_string(directory_sexp, "directory");

    f->read_from_directory(directory);
  }

  void Forwarder_read_from_directory_and_no_states(SEXP f_ptr_sexp, SEXP directory_sexp, SEXP no_states_sexp) {
    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);
    std::string directory = parse_string(directory_sexp, "directory");
    size_t no_states = parse_integer(no_states_sexp, "no_states");

    f->read_from_directory(directory, no_states);
  }

  void Forwarder_write_to_directory(SEXP f_ptr_sexp, SEXP directory_sexp) {
    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);
    std::string directory = parse_string(directory_sexp, "directory");

    f->write_to_directory(directory);
  }

  SEXP Forwarder_forward(SEXP f_ptr_sexp, SEXP pi_sexp, SEXP A_sexp, SEXP B_sexp, SEXP no_states_sexp, SEXP alphabet_size_sexp) { 
    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);

    zipHMM::Matrix pi, A, B;
    convert_HMM_parameters(pi_sexp, A_sexp, B_sexp, no_states_sexp, alphabet_size_sexp, pi, A, B);

    double loglik = f->forward(pi, A, B);

    SEXP loglik_sexp;
    PROTECT(loglik_sexp = NEW_NUMERIC(1.0));
    NUMERIC_POINTER(loglik_sexp)[0] = loglik;

    UNPROTECT(1);
    return loglik_sexp;
  }

  SEXP Forwarder_ptforward(SEXP f_ptr_sexp, SEXP pi_sexp, SEXP A_sexp, SEXP B_sexp, SEXP no_states_sexp, SEXP alphabet_size_sexp) { 

    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);

    zipHMM::Matrix pi, A, B;
    convert_HMM_parameters(pi_sexp, A_sexp, B_sexp, no_states_sexp, alphabet_size_sexp, pi, A, B);

    double loglik = f->pthread_forward(pi, A, B);

    SEXP loglik_sexp;
    PROTECT(loglik_sexp = NEW_NUMERIC(1.0));
    NUMERIC_POINTER(loglik_sexp)[0] = loglik;

    UNPROTECT(1);
    return loglik_sexp;
  }

  SEXP Forwarder_get_orig_seq_length(SEXP f_ptr_sexp) { 

    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);

    size_t orig_seq_length = f->get_orig_seq_length();

    SEXP orig_seq_length_sexp;
    PROTECT(orig_seq_length_sexp = NEW_INTEGER(1));
    INTEGER_POINTER(orig_seq_length_sexp)[0] = orig_seq_length;

    UNPROTECT(1);
    return orig_seq_length_sexp;
  }

  SEXP Forwarder_get_orig_alphabet_size(SEXP f_ptr_sexp) { 

    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);

    size_t orig_alphabet_size = f->get_orig_alphabet_size();

    SEXP orig_alphabet_size_sexp;
    PROTECT(orig_alphabet_size_sexp = NEW_INTEGER(1));
    INTEGER_POINTER(orig_alphabet_size_sexp)[0] = orig_alphabet_size;

    UNPROTECT(1);
    return orig_alphabet_size_sexp;
  }

  SEXP Forwarder_get_seq_length(SEXP f_ptr_sexp, SEXP no_states_sexp) { 

    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);

    size_t no_states = parse_integer(no_states_sexp, "no_states");

    size_t seq_length = f->get_seq_length(no_states);

    SEXP seq_length_sexp;
    PROTECT(seq_length_sexp = NEW_INTEGER(1));
    INTEGER_POINTER(seq_length_sexp)[0] = seq_length;

    UNPROTECT(1);
    return seq_length_sexp;
  }

  SEXP Forwarder_get_alphabet_size(SEXP f_ptr_sexp, SEXP no_states_sexp) { 

    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);

    size_t no_states = parse_integer(no_states_sexp, "no_states");

    size_t alphabet_size = f->get_alphabet_size(no_states);

    SEXP alphabet_size_sexp;
    PROTECT(alphabet_size_sexp = NEW_INTEGER(1));
    INTEGER_POINTER(alphabet_size_sexp)[0] = alphabet_size;

    UNPROTECT(1);
    return alphabet_size_sexp;
  }

  SEXP Forwarder_get_pair(SEXP f_ptr_sexp, SEXP symbol_sexp) {

    zipHMM::Forwarder *f = cast_forwarder(f_ptr_sexp);

    size_t symbol = parse_integer(symbol_sexp, "symbol");

    zipHMM::s_pair pair = f->get_pair(symbol); 
    
    SEXP pair_sexp = PROTECT(allocVector(INTSXP, 2));
    INTEGER(pair_sexp)[0] = pair.first;
    INTEGER(pair_sexp)[1] = pair.second;

    UNPROTECT(1);
    return pair_sexp;
  }

}

// Viterbi
extern "C" {
  SEXP c_viterbi(SEXP seq_filename_sexp, SEXP pi_sexp, SEXP A_sexp, SEXP B_sexp, SEXP no_states_sexp, SEXP alphabet_size_sexp) {
    std::string seq_filename = parse_string(seq_filename_sexp, "seq_filename");

    double *pi_content, *A_content, *B_content;
    zipHMM::Matrix c_pi, c_A, c_B;
    parse_HMM_parameters(pi_sexp, A_sexp, B_sexp, &pi_content, &A_content, &B_content);
    convert_HMM_parameters(pi_sexp, A_sexp, B_sexp, no_states_sexp, alphabet_size_sexp, c_pi, c_A, c_B);
    
    std::vector<unsigned> viterbi_path;
    double ll = viterbi(seq_filename, c_pi, c_A, c_B, viterbi_path);

    SEXP ll_sexp;
    PROTECT(ll_sexp = NEW_NUMERIC(1.0));
    NUMERIC_POINTER(ll_sexp)[0] = ll;

    // parse path
    SEXP viterbi_path_sexp;
    PROTECT(viterbi_path_sexp = NEW_INTEGER(viterbi_path.size()));
    for(size_t i = 0; i < viterbi_path.size(); ++i)
      INTEGER(viterbi_path_sexp)[i] = viterbi_path[i];
    
    SEXP out;
    PROTECT(out = NEW_LIST(2));
    SET_VECTOR_ELT(out, 0, ll_sexp);
    SET_VECTOR_ELT(out, 1, viterbi_path_sexp);

    UNPROTECT(3);
    return out;
  }
}

// Posterior decoding
extern "C" {
  SEXP c_posterior_decoding(SEXP seq_filename_sexp, SEXP pi_sexp, SEXP A_sexp, SEXP B_sexp, SEXP no_states_sexp, SEXP alphabet_size_sexp) {
    std::string seq_filename = parse_string(seq_filename_sexp, "seq_filename");

    double *pi_content, *A_content, *B_content;
    zipHMM::Matrix c_pi, c_A, c_B;
    parse_HMM_parameters(pi_sexp, A_sexp, B_sexp, &pi_content, &A_content, &B_content);
    convert_HMM_parameters(pi_sexp, A_sexp, B_sexp, no_states_sexp, alphabet_size_sexp, c_pi, c_A, c_B);

    std::vector<unsigned> pd_path;
    zipHMM::Matrix pd_table;
    posterior_decoding(seq_filename, c_pi, c_A, c_B, pd_path, pd_table);

    SEXP pd_path_sexp;
    PROTECT(pd_path_sexp = NEW_INTEGER(pd_path.size()));
    for(size_t i = 0; i < pd_path.size(); ++i)
      INTEGER(pd_path_sexp)[i] = pd_path[i];

    SEXP pd_table_sexp = PROTECT(allocMatrix(REALSXP, pd_table.get_height(), pd_table.get_width()));
    double *pd_table_content = REAL(pd_table_sexp);
    for(size_t r = 0; r < pd_table.get_height(); ++r)
      for(size_t c = 0; c < pd_table.get_width(); ++c)
	pd_table_content[c*pd_table.get_height() + r] = pd_table(r, c);

    SEXP out;
    PROTECT(out = NEW_LIST(2));
    
    SET_VECTOR_ELT(out, 0, pd_path_sexp);
    SET_VECTOR_ELT(out, 1, pd_table_sexp);

    UNPROTECT(3);
    return out;
  }
}

// calibrate
extern "C" {
  void c_calibrate(SEXP device_filename_sexp) {
    std::string device_filename = parse_string(device_filename_sexp, "device_filename");
    if(std::strcmp(device_filename.c_str(), "-") == 0)
      device_filename = zipHMM::DEFAULT_DEVICE_FILENAME;

    zipHMM::calibrate(device_filename);
  }
}
