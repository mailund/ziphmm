#include "forwarder.hpp"
#include "simple_forwarder.hpp"
#include "matrix.hpp"
#include "hmm_io.hpp"
#include "posterior_decoding.hpp"
#include "viterbi.hpp"
#include "calibrate.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <new>
#include <Python.h>

// HMM_IO
extern "C" {
  int c_read_HMM_spec(size_t &no_states, size_t &alphabet_size, const char *filename) {
    zipHMM::read_HMM_spec(no_states, alphabet_size, std::string(filename));
    return 1;
  }

  int c_read_HMM(void *pi_ptr, void *A_ptr, void *B_ptr, const char *filename) {
    zipHMM::Matrix *pi = reinterpret_cast<zipHMM::Matrix *>(pi_ptr);
    zipHMM::Matrix *A = reinterpret_cast<zipHMM::Matrix *>(A_ptr);
    zipHMM::Matrix *B = reinterpret_cast<zipHMM::Matrix *>(B_ptr);

    zipHMM::read_HMM(*pi, *A, *B, std::string(filename));
    return 1;
  }

  int c_write_HMM(void *pi_ptr, void *A_ptr, void *B_ptr, const char *filename) {
    zipHMM::Matrix *pi = reinterpret_cast<zipHMM::Matrix *>(pi_ptr);
    zipHMM::Matrix *A = reinterpret_cast<zipHMM::Matrix *>(A_ptr);
    zipHMM::Matrix *B = reinterpret_cast<zipHMM::Matrix *>(B_ptr);
    
    zipHMM::write_HMM(*pi, *A, *B, std::string(filename));
    return 1;
  }

  PyObject *c_read_seq(const char *seq_filename) {
    std::vector<unsigned> sequence;
    zipHMM::readSeq(sequence, seq_filename);

    Py_Initialize();

    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();
    
    PyObject *result = PyList_New(0);
    for(std::vector<unsigned>::const_iterator it = sequence.begin(); it != sequence.end(); ++it) {
      PyObject *i = PyLong_FromLong(*it);
      PyList_Append(result, i);
      Py_DECREF(i);
    }
    
    PyGILState_Release(gstate);  

    return result;
  }


}

// Forwarder
extern "C" {
  
  void *Forwarder_new() {
    return new zipHMM::Forwarder();
  }

  void Forwarder_read_seq(void *f_ptr, const char *seq_filename, const size_t alphabet_size, unsigned *nStatesSave_array, size_t nStatesSave_length, size_t min_no_evals) {
    zipHMM::Forwarder *forwarder = reinterpret_cast<zipHMM::Forwarder *>(f_ptr);

    std::vector<size_t> nStatesSave_vector;
    for(size_t i = 0; i < nStatesSave_length; ++i)
      nStatesSave_vector.push_back(nStatesSave_array[i]);

    forwarder->read_seq(seq_filename, alphabet_size, nStatesSave_vector, min_no_evals); 
  }

  // void read_seq_directory(const std::string &dir_filename, const size_t alphabet_size, std::vector<size_t> &nStatesSave, const size_t min_no_eval = 1);

  void Forwarder_read_seq_directory(void *f_ptr, const char *dirname, const size_t alphabet_size, unsigned *nStatesSave_array, size_t nStatesSave_length, const size_t min_no_evals) {
    zipHMM::Forwarder *forwarder = reinterpret_cast<zipHMM::Forwarder *>(f_ptr);

    std::vector<size_t> nStatesSave_vector;
    for(size_t i = 0; i < nStatesSave_length; ++i)
      nStatesSave_vector.push_back(nStatesSave_array[i]);

    forwarder->read_seq_directory(dirname, alphabet_size, nStatesSave_vector, min_no_evals);
  }

  void Forwarder_read_from_directory(void *f_ptr, const char *directory) {
    zipHMM::Forwarder *forwarder = reinterpret_cast<zipHMM::Forwarder *>(f_ptr);
    forwarder->read_from_directory(directory);
  }

  void Forwarder_read_from_directory_and_no_states(void *f_ptr, const char *directory, const size_t no_states) {
    zipHMM::Forwarder *forwarder = reinterpret_cast<zipHMM::Forwarder *>(f_ptr);
    forwarder->read_from_directory(directory, no_states);
  }

  int Forwarder_destructor(void *f_ptr) {
    zipHMM::Forwarder *f = reinterpret_cast<zipHMM::Forwarder *>(f_ptr);
    delete f;
    return 1;
  }

  double Forwarder_forward(const void *f_ptr, const void *pi_ptr, const void *A_ptr, const void *B_ptr) { 
    const zipHMM::Forwarder *f = reinterpret_cast<const zipHMM::Forwarder *>(f_ptr);
    const zipHMM::Matrix *pi = reinterpret_cast<const zipHMM::Matrix *>(pi_ptr);
    const zipHMM::Matrix *A = reinterpret_cast<const zipHMM::Matrix *>(A_ptr);
    const zipHMM::Matrix *B = reinterpret_cast<const zipHMM::Matrix *>(B_ptr);
    return f->forward(*pi, *A, *B); 
  }

  double Forwarder_pthread_forward(const void *f_ptr, const void *pi_ptr, const void *A_ptr, const void *B_ptr, const char *device_filename) { 
    const zipHMM::Forwarder *f = reinterpret_cast<const zipHMM::Forwarder *>(f_ptr);
    const zipHMM::Matrix *pi = reinterpret_cast<const zipHMM::Matrix *>(pi_ptr);
    const zipHMM::Matrix *A = reinterpret_cast<const zipHMM::Matrix *>(A_ptr);
    const zipHMM::Matrix *B = reinterpret_cast<const zipHMM::Matrix *>(B_ptr);
    if(strcmp(device_filename, "-") == 0)
      return f->pthread_forward(*pi, *A, *B);
    else
      return f->pthread_forward(*pi, *A, *B, device_filename);
  }

  double Forwarder_mr_pthread_forward(const void *f_ptr, const void *pi_ptr, const void *A_ptr, const void *B_ptr, const char *device_filename) { 
    const zipHMM::Forwarder *f = reinterpret_cast<const zipHMM::Forwarder *>(f_ptr);
    const zipHMM::Matrix *pi = reinterpret_cast<const zipHMM::Matrix *>(pi_ptr);
    const zipHMM::Matrix *A = reinterpret_cast<const zipHMM::Matrix *>(A_ptr);
    const zipHMM::Matrix *B = reinterpret_cast<const zipHMM::Matrix *>(B_ptr);
    if(strcmp(device_filename, "-") == 0)
      return f->mr_pthread_forward(*pi, *A, *B);
    else
      return f->mr_pthread_forward(*pi, *A, *B, device_filename);
  }

  size_t Forwarder_get_orig_seq_length(const void *f_ptr) { 
    const zipHMM::Forwarder *f = reinterpret_cast<const zipHMM::Forwarder *>(f_ptr);
    return f->get_orig_seq_length();
  }

  size_t Forwarder_get_orig_alphabet_size(const void *f_ptr) { 
    const zipHMM::Forwarder *f = reinterpret_cast<const zipHMM::Forwarder *>(f_ptr);
    return f->get_orig_alphabet_size(); 
  }
  
  size_t Forwarder_get_seq_length(const void *f_ptr, size_t no_states) { 
    const zipHMM::Forwarder *f = reinterpret_cast<const zipHMM::Forwarder *>(f_ptr);
    return f->get_seq_length(no_states); 
  }

  size_t Forwarder_get_alphabet_size(const void *f_ptr, size_t no_states) { 
    const zipHMM::Forwarder *f = reinterpret_cast<const zipHMM::Forwarder *>(f_ptr);
    return f->get_alphabet_size(no_states); 
  }

  PyObject *Forwarder_get_pair(const void *f_ptr, unsigned symbol) {
    const zipHMM::Forwarder *f = reinterpret_cast<const zipHMM::Forwarder *>(f_ptr);
    zipHMM::s_pair pair = f->get_pair(symbol); 
    return Py_BuildValue("ll", pair.first, pair.second);
  }

  void Forwarder_write_to_directory(const void *f_ptr, const char *directory) {
    const zipHMM::Forwarder *f = reinterpret_cast<const zipHMM::Forwarder *>(f_ptr);
    f->write_to_directory(directory);
  }
}

// SimpleForwarder

extern "C" {

  void *SimpleForwarder_new(const char *seq_filename) {
    return new zipHMM::SimpleForwarder(seq_filename);
  }

  double SimpleForwarder_forward(const void *f_ptr, const void *pi_ptr, const void *A_ptr, const void *B_ptr) { 
    const zipHMM::SimpleForwarder *f = reinterpret_cast<const zipHMM::SimpleForwarder *>(f_ptr);
    const zipHMM::Matrix *pi = reinterpret_cast<const zipHMM::Matrix *>(pi_ptr);
    const zipHMM::Matrix *A = reinterpret_cast<const zipHMM::Matrix *>(A_ptr);
    const zipHMM::Matrix *B = reinterpret_cast<const zipHMM::Matrix *>(B_ptr);
    return f->forward(*pi, *A, *B); 
  }

  int SimpleForwarder_destructor(void *f_ptr) {
    zipHMM::SimpleForwarder *f = reinterpret_cast<zipHMM::SimpleForwarder *>(f_ptr);
    delete f;
    return 1;
  }

}

// posterior decoding
extern "C" {
  int c_posterior_decoding(void *pd_path_vptr, void *pd_table_vptr,
			   void *pi_vptr, void *A_vptr, void *B_vptr,
			   const char *seq_filename) {
    std::vector<unsigned> *pd_path_ptr = reinterpret_cast<std::vector<unsigned> *>(pd_path_vptr);
    zipHMM::Matrix *pd_table_ptr = reinterpret_cast<zipHMM::Matrix *>(pd_table_vptr);
    zipHMM::Matrix *pi_ptr = reinterpret_cast<zipHMM::Matrix *>(pi_vptr);
    zipHMM::Matrix *A_ptr = reinterpret_cast<zipHMM::Matrix *>(A_vptr);
    zipHMM::Matrix *B_ptr = reinterpret_cast<zipHMM::Matrix *>(B_vptr);
    
    zipHMM::posterior_decoding(seq_filename, *pi_ptr, *A_ptr, *B_ptr, *pd_path_ptr, *pd_table_ptr);
    
    return 1;
  }
}

// Viterbi
extern "C" {
  double c_viterbi(void *viterbi_path_vptr,
		   void *pi_vptr, void *A_vptr, void *B_vptr,
		   const char *seq_filename) {
    std::vector<unsigned> *viterbi_path_ptr = reinterpret_cast<std::vector<unsigned> *>(viterbi_path_vptr);
    zipHMM::Matrix *pi_ptr = reinterpret_cast<zipHMM::Matrix *>(pi_vptr);
    zipHMM::Matrix *A_ptr = reinterpret_cast<zipHMM::Matrix *>(A_vptr);
    zipHMM::Matrix *B_ptr = reinterpret_cast<zipHMM::Matrix *>(B_vptr);
    
    double viterbi_ll = zipHMM::viterbi(seq_filename, *pi_ptr, *A_ptr, *B_ptr, *viterbi_path_ptr);
    
    return viterbi_ll;
  }
}



// Sequence
extern "C" {
  void *Sequence_new() {
    return new(std::nothrow) std::vector<unsigned>();
  }

  int Sequence_destructor(void *seq_ptr) {
    std::vector<unsigned> *seq = reinterpret_cast< std::vector<unsigned> * >(seq_ptr);
    delete seq;
    return 1;
  }

  unsigned Sequence_get(const void *seq_ptr, size_t i) { 
    const std::vector<unsigned> *seq = reinterpret_cast<const std::vector<unsigned> *>(seq_ptr);
    return (*seq)[i];
  }

  unsigned Sequence_len(const void *seq_ptr) {
    const std::vector<unsigned> *seq = reinterpret_cast<const std::vector<unsigned> *>(seq_ptr);
    return unsigned(seq->size());
  }
}

// Matrix
extern "C" {
  void *Matrix_new_empty() { 
    return new(std::nothrow) zipHMM::Matrix(); 
  }

  void *Matrix_new_height_width(size_t height, size_t width) { 
    return new(std::nothrow) zipHMM::Matrix(height, width); 
  }

  int Matrix_destructor(void *m_ptr) {
    zipHMM::Matrix *m = reinterpret_cast<zipHMM::Matrix *>(m_ptr);
    delete m;
    return 1;
  }

  size_t Matrix_get_width(const void *m_ptr) { 
    const zipHMM::Matrix *m = reinterpret_cast<const zipHMM::Matrix *>(m_ptr);
    return m->get_width(); 
  }

  size_t Matrix_get_height(const void *m_ptr) { 
    const zipHMM::Matrix *m = reinterpret_cast<const zipHMM::Matrix *>(m_ptr);
    return m->get_height();
  }

  int Matrix_reset(void *m_ptr, size_t height, size_t width) { 
    zipHMM::Matrix *m = reinterpret_cast<zipHMM::Matrix *>(m_ptr);
    m->reset(height, width); 
    return 1;
  }
  
  int Matrix_set(void *m_ptr, size_t i, size_t j, double value) { 
    zipHMM::Matrix *m = reinterpret_cast<zipHMM::Matrix *>(m_ptr);
    m->set(i, j, value); 
    return 1;
  }

  double Matrix_get(const void *m_ptr, size_t i, size_t j) { 
    const zipHMM::Matrix *m = reinterpret_cast<const zipHMM::Matrix *>(m_ptr);
    return m->get(i, j); 
  }

  int Matrix_transpose(const void *from_ptr, void *to_ptr) {
    const zipHMM::Matrix *from = reinterpret_cast<const zipHMM::Matrix *>(from_ptr);
    zipHMM::Matrix *to = reinterpret_cast<zipHMM::Matrix *>(to_ptr);

    zipHMM::Matrix::transpose(*from, *to);
    return 1;
  }

  int Matrix_print(const void *m_ptr) { 
    const zipHMM::Matrix *m = reinterpret_cast<const zipHMM::Matrix *>(m_ptr);
    m->print();
    return 1;
  }
}

// Calibrate
extern "C" {
  int c_calibrate(const char *device_filename) {
    zipHMM::calibrate(device_filename);

    return 0;
  }
}
