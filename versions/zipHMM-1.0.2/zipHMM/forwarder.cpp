#include "forwarder.hpp"
#include "matrix.hpp"
#include "seq_io.hpp"
#include "io_utils.hpp"
#include "hmm_utils.hpp"
#include "PThreadProcessingDevice.hpp"
#include "performance_description.hpp"
#include "Stage1JobControl.hpp"
#include "Stage2JobControl.hpp"
#include "debug.hpp"

#include <utility>
#include <algorithm>
#include <map>
#include <vector>
#include <deque>
#include <stack>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <climits>
#include <dirent.h>

#ifdef WITH_CUDA
#include "cutil.h"
#include <cuda_runtime.h> 
#include "cublas_v2.h"
#endif


namespace zipHMM {

  void Forwarder::read_seq(const std::string &seq_filename, const size_t alphabet_size, std::vector<size_t> nStatesSave, const size_t min_no_eval) {
    orig_alphabet_size = alphabet_size;
    std::vector<unsigned> *prev_seq_p = 0;
    std::vector<unsigned> *current_seq_p = 0;
    std::vector<unsigned> *next_seq_p = 0;
    std::vector<size_t> pair_n;
    bool counted = false;
    size_t prev_max_count = UINT_MAX;
    Timer iteration_timer;
    bool nStatesSave_init_empty = nStatesSave.empty();
    size_t no_states_save = 0;
    std::map<size_t, double> nStates2t_MM;
    std::map<size_t, double> nStates2t_MV;
    double iteration_time;

    if(!nStatesSave_init_empty) {
      no_states_save = nStatesSave.back();
      for(std::vector<size_t>::iterator it = nStatesSave.begin(); it != nStatesSave.end(); ++it) {
	nStates2t_MM[*it] = Matrix::time_blas_mult((*it));
	nStates2t_MV[*it] = Matrix::time_blas_matrix_vector_mult((*it));
      }
    }
    
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

    iteration_timer.start();

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

      iteration_timer.stop();
      iteration_time = iteration_timer.timeElapsed();

      // stopping criteria
      if(nStatesSave_init_empty && prev_max_count == max_count) {
	nStates2seq[2] = *current_seq_p; // copying sequence
	nStates2alphabet_size[2] = new_alphabet_size;
	break;
      } else if(!nStatesSave_init_empty) {
	double time_saved_per_eval = double(max_count) * double(nStates2t_MV[no_states_save]) - double(nStates2t_MM[no_states_save]);
	double time_saved_on_evals = double(min_no_eval) * time_saved_per_eval;
	double total_time_saved = time_saved_on_evals - iteration_time;
	
	if(total_time_saved <= 0) {
	  nStates2seq[no_states_save] = *current_seq_p; // copying sequence
	  nStates2alphabet_size[no_states_save] = new_alphabet_size;

	  nStatesSave.pop_back();
	  if(nStatesSave.empty())
	    break;
	  else
	    no_states_save = nStatesSave.back();
	}
      }
      
      iteration_timer.start();

      // save the components of max_pair in symbol2pair
      symbol2pair[unsigned(new_alphabet_size)] = max_pair;

      pair_n.clear();
      pair_n.resize( (new_alphabet_size+1) * (new_alphabet_size+1));
      next_seq_p = new std::vector<unsigned>();
      next_seq_p->push_back( (*current_seq_p)[0] );
      // replace every occurrence of max_pair with a new alphabet symbol and count pairs at the same time
      size_t i = 1; // first position is special because we have not seen any pairs yet.
      if(matches(max_pair, current_seq_p, i)) { 
	next_seq_p->push_back(unsigned(new_alphabet_size));
	++i;
      } else {
	next_seq_p->push_back((*current_seq_p)[i]);
      }
      // do the rest of the sequence
      for(i = i+1; i < current_seq_p->size() - 1; i++) {
	if(matches(max_pair, current_seq_p, i)) {
	  next_seq_p->push_back(unsigned(new_alphabet_size));
	  ++i;
	} else {
	  next_seq_p->push_back((*current_seq_p)[i]);
	}

	add_count(pair_n, next_seq_p, next_seq_p->size() - 2, counted, new_alphabet_size + 1);
      }
      if(i == current_seq_p->size() - 1) { // push last symbol of sequence
	next_seq_p->push_back((*current_seq_p)[i]);
	add_count(pair_n, next_seq_p, next_seq_p->size() - 2, counted, new_alphabet_size + 1);
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

  void Forwarder::read_seq(const std::string &seq_filename, const size_t alphabet_size, const size_t no_states, const size_t min_no_eval) {
    std::vector<size_t> nStatesSave;
    nStatesSave.push_back(no_states);
    read_seq(seq_filename, alphabet_size, nStatesSave, min_no_eval);
  }

  void Forwarder::read_seq(const std::string &seq_filename, const size_t alphabet_size, const size_t min_no_eval) {
    std::vector<size_t> nStatesSave;
    read_seq(seq_filename, alphabet_size, nStatesSave, min_no_eval);
  }

  void Forwarder::read_from_directory(const std::string &directory) {
    read_data_structure_from_directory(directory);
    read_seqs_from_directory(directory + "/nStates2seq");
  }

  void Forwarder::read_from_directory(const std::string &directory, const size_t no_states) {
    read_data_structure_from_directory(directory);
    std::map<size_t, size_t> new_nStates2alphabet_size;
    // remove what we read too much
    new_nStates2alphabet_size[no_states] = nStates2alphabet_size[no_states];
    nStates2alphabet_size = new_nStates2alphabet_size;

    read_seq_from_directory(directory + "/nStates2seq", no_states);    
  }

  double Forwarder::forward(const Matrix &pi, const Matrix &A, const Matrix &B) const {
    if(pi.get_width() != 1 || pi.get_height() != A.get_width() || A.get_height() != A.get_width() ||
       B.get_height() != A.get_width() || B.get_width() != orig_alphabet_size) {
      std::cerr << "Dimensions of input matrices do not match:" << std::endl;
      std::cerr << "\t" << "pi width:  " << pi.get_width()  << std::endl;
      std::cerr << "\t" << "pi height: " << pi.get_height() << std::endl;
      std::cerr << "\t" << "A width:  "  << A.get_width()   << std::endl;
      std::cerr << "\t" << "A height: "  << A.get_height()  << std::endl;
      std::cerr << "\t" << "B width:  "  << B.get_width()   << std::endl;
      std::cerr << "\t" << "B height: "  << B.get_height()  << std::endl;
      std::exit(-1);
    }

    const std::vector<unsigned> *seq = 0;
    std::deque<double> scales;
    double *symbol2scale;
    Matrix *symbol2matrix;
    Matrix res;
    Matrix tmp;
    double loglikelihood = 0;

    // find alphabet and seq for given number of states
    size_t no_states = A.get_width();
    size_t alphabet_size = 0;
    for(std::map<size_t, std::vector<unsigned> >::const_iterator it = nStates2seq.begin(); it != nStates2seq.end(); ++it) {
      if(it->first >= no_states) {
	seq = &(it->second);
	alphabet_size = nStates2alphabet_size.find(it->first)->second;
	break;
      }
    }
    if(seq == 0) {
      seq = &(nStates2seq.rbegin()->second);
      alphabet_size = nStates2alphabet_size.rbegin()->second;
    }
    
    symbol2scale = new double[alphabet_size];
    symbol2matrix = new Matrix[alphabet_size];    

    compute_symbol2scale_and_symbol2matrix(symbol2matrix, symbol2scale, A, B, alphabet_size);

    // compute C_1 and push corresponding scale
    scales.push_back( std::log(init_apply_em_prob(res, pi, B, (*seq)[0])) );

    // multiply matrices across the sequence
    for(size_t i = 1; i < seq->size(); ++i) {
      Matrix::blas_matrix_vector_mult(symbol2matrix[(*seq)[i]], res, tmp);
      Matrix::copy(tmp, res);

      scales.push_back( std::log(res.normalize()) );
      scales.push_back( symbol2scale[(*seq)[i]] );
    }
   
    // compute loglikelihood by summing log of scales
    for(std::deque<double>::iterator it = scales.begin(); it != scales.end(); ++it) {
      loglikelihood += (*it);
    }

    delete[] symbol2scale;
    delete[] symbol2matrix;

    return loglikelihood;
  }

  double Forwarder::pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B, const DeviceDescriptor &device_descriptor) const {
    if(pi.get_width() != 1 || pi.get_height() != A.get_width() || A.get_height() != A.get_width() ||
       B.get_height() != A.get_width() || B.get_width() != orig_alphabet_size) {
      std::cerr << "Dimensions of input matrices do not match:" << std::endl;
      std::cerr << "\t" << "pi width:  " << pi.get_width()  << std::endl;
      std::cerr << "\t" << "pi height: " << pi.get_height() << std::endl;
      std::cerr << "\t" << "A width:  "  << A.get_width()   << std::endl;
      std::cerr << "\t" << "A height: "  << A.get_height()  << std::endl;
      std::cerr << "\t" << "B width:  "  << B.get_width()   << std::endl;
      std::cerr << "\t" << "B height: "  << B.get_height()  << std::endl;
      std::exit(-1);
    }

    const std::vector<unsigned> *seq = 0;
    std::deque<double> scales;
    double *symbol2scale;
    Matrix *symbol2matrix;
    std::vector<ProcessingDevice*> devices;
    size_t numberOfDevices;
    Matrix *result;
    Matrix *temp;
    double loglikelihood;
    size_t head;

    // find alphabet and seq for given number of states
    size_t no_states = A.get_width();
    size_t alphabet_size = 0;
    for(std::map<size_t, std::vector<unsigned> >::const_iterator it = nStates2seq.begin(); it != nStates2seq.end(); ++it) {
      if(it->first >= no_states) {
	seq = &(it->second);
	alphabet_size = nStates2alphabet_size.find(it->first)->second;
	break;
      }
    }
    if(seq == 0) {
      seq = &(nStates2seq.rbegin()->second);
      alphabet_size = nStates2alphabet_size.rbegin()->second;
    }

    const size_t length = seq->size();
    const size_t nBlocks = size_t(std::sqrt(length));
    Stage2JobControl control(length, nBlocks);
    symbol2scale = new double[alphabet_size];
    symbol2matrix = new Matrix[alphabet_size];

    compute_symbol2scale_and_symbol2matrix(symbol2matrix, symbol2scale, A, B, alphabet_size);

    numberOfDevices = device_descriptor.getNDevices();
    
    for(unsigned i = 0; i < numberOfDevices; ++i) {
      devices.push_back(device_descriptor.createDevice(i));
      devices[i]->setParameters(&pi, &A, &B, &symbol2pair, symbol2scale, symbol2matrix, seq);
    }
    
    devices[0]->likelihoodVector(control);
    for(unsigned i = 1; i < numberOfDevices; ++i)
      devices[i]->likelihoodMatrix(control);
    
    for(unsigned i = 0; i < numberOfDevices; ++i)
      devices[i]->join();

    result = new Matrix();
    temp = new Matrix();
  
    head = control.headBlock - 1;
  
    Matrix::copy(*control.resultMatrices[head], *result);
    loglikelihood = control.resultLogLikelihoods[head];
  
    for(size_t i = head + 1; i < nBlocks; ++i) {
      Matrix::mult(*control.resultMatrices[i], *result, *temp);
      std::swap(result, temp);
    
      loglikelihood += LinearSpace::toLogSpace(result->normalize()) + control.resultLogLikelihoods[i];
    }

    for(unsigned i = 0; i < devices.size(); ++i)
      delete devices[i];
  
    delete result;
    delete temp;
  
    delete[] symbol2scale;
    delete[] symbol2matrix;

    return loglikelihood;
  }

  double Forwarder::pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B, const std::string &device_filename) const {
    std::vector<DeviceDescriptor> device_descriptors;
    readDescriptors(device_descriptors, device_filename);
    return pthread_forward(pi, A, B, device_descriptors[0]);
  }

  double Forwarder::pthread_forward_par_stage1(const Matrix &pi, const Matrix &A, const Matrix &B, const DeviceDescriptor &device_descriptor) const {
    if(pi.get_width() != 1 || pi.get_height() != A.get_width() || A.get_height() != A.get_width() ||
       B.get_height() != A.get_width() || B.get_width() != orig_alphabet_size) {
      std::cerr << "Dimensions of input matrices do not match:" << std::endl;
      std::cerr << "\t" << "pi width:  " << pi.get_width()  << std::endl;
      std::cerr << "\t" << "pi height: " << pi.get_height() << std::endl;
      std::cerr << "\t" << "A width:  "  << A.get_width()   << std::endl;
      std::cerr << "\t" << "A height: "  << A.get_height()  << std::endl;
      std::cerr << "\t" << "B width:  "  << B.get_width()   << std::endl;
      std::cerr << "\t" << "B height: "  << B.get_height()  << std::endl;
      std::exit(-1);
    }

    const std::vector<unsigned> *seq = 0;
    std::deque<double> scales;
    std::vector<ProcessingDevice*> devices;
    size_t numberOfDevices;
    Matrix *result;
    Matrix *temp;
    double loglikelihood;
    size_t head;

    // find alphabet and seq for given number of states
    size_t no_states = A.get_width();
    size_t alphabet_size = 0;
    for(std::map<size_t, std::vector<unsigned> >::const_iterator it = nStates2seq.begin(); it != nStates2seq.end(); ++it) {
      if(it->first >= no_states) {
	seq = &(it->second);
	alphabet_size = nStates2alphabet_size.find(it->first)->second;
	break;
      }
    }
    if(seq == 0) {
      seq = &(nStates2seq.rbegin()->second);
      alphabet_size = nStates2alphabet_size.rbegin()->second;
    }

    const size_t length = seq->size();
    const size_t nBlocks = size_t(std::sqrt(length));
    Stage2JobControl control(length, nBlocks);
    double *symbol2scale = new double[alphabet_size];
    Matrix *symbol2matrix = new Matrix[alphabet_size];

    numberOfDevices = device_descriptor.getNDevices();

    for(unsigned i = 0; i < numberOfDevices; ++i) {
      devices.push_back(device_descriptor.createDevice(i));
      devices[i]->setParameters(&pi, &A, &B, &symbol2pair, symbol2scale, symbol2matrix, seq);
    }

    // stage 1
    pthread_compute_symbol2scale_and_symbol2matrix(symbol2scale, symbol2matrix, A, B, alphabet_size, devices);

    // stage 2    
    devices[0]->likelihoodVector(control);
    for(unsigned i = 1; i < numberOfDevices; ++i)
      devices[i]->likelihoodMatrix(control);
    
    for(unsigned i = 0; i < numberOfDevices; ++i)
      devices[i]->join();

    // stage 3
    result = new Matrix();
    temp = new Matrix();
  
    head = control.headBlock - 1;
  
    Matrix::copy(*control.resultMatrices[head], *result);
    loglikelihood = control.resultLogLikelihoods[head];
  
    for(size_t i = head + 1; i < nBlocks; ++i) {
      Matrix::mult(*control.resultMatrices[i], *result, *temp);
      std::swap(result, temp);
    
      loglikelihood += LinearSpace::toLogSpace(result->normalize()) + control.resultLogLikelihoods[i];
    }

    // clean up
    for(unsigned i = 0; i < devices.size(); ++i)
      delete devices[i];
  
    delete result;
    delete temp;
  
    delete[] symbol2scale;
    delete[] symbol2matrix;

    return loglikelihood;
  }

  double Forwarder::pthread_forward_par_stage1(const Matrix &pi, const Matrix &A, const Matrix &B, const std::string &device_filename) const {
    std::vector<DeviceDescriptor> device_descriptors;
    readDescriptors(device_descriptors, device_filename);
    return pthread_forward_par_stage1(pi, A, B, device_descriptors[0]);
  }


  void Forwarder::write_to_directory(const std::string &directory) const {
    std::string wd = get_working_directory();
    std::string absolute_dir_name = wd + "/" + directory;
    std::string data_structure_filename = absolute_dir_name + "/data_structure";
    std::string nStates2seq_absolute_dir_name = absolute_dir_name + "/nStates2seq";
    
    // create directory
    mk_dir(absolute_dir_name);
    mk_dir(nStates2seq_absolute_dir_name);

    // write basic data structure
    write_data_structure(data_structure_filename);

    // write nStates2seq
    write_seqs(nStates2seq_absolute_dir_name);
  }

  void Forwarder::write_seqs(const std::string nstates2seq_absolute_dir_name) const {
    for(std::map<size_t, std::vector<unsigned> >::const_iterator it = nStates2seq.begin(); it != nStates2seq.end(); ++it) {
      size_t nStates = (*it).first;
      
      std::stringstream seq_filename_stream;
      seq_filename_stream << nstates2seq_absolute_dir_name << "/" << nStates << ".seq";
      
      std::string seq_filename = seq_filename_stream.str();
      write_seq(seq_filename, nStates);
    }
  }

  void Forwarder::write_seq(const std::string &seq_filename, size_t no_states) const {
    std::ofstream seq_stream(seq_filename.c_str());
    if(!seq_stream) {
      std::cerr << "Unable to open \"" << seq_filename << "\"" << std::endl;
      exit(-1);
    }

    write_seq(seq_stream, no_states);
  }

  void Forwarder::write_seq(std::ofstream &out, size_t no_states) const {
    const std::vector<unsigned> &nStatesSeq = nStates2seq.find(no_states)->second;
    for(std::vector<unsigned>::const_iterator it = nStatesSeq.begin(); it != nStatesSeq.end(); ++it)
      out << (*it) << " ";
  }

  void Forwarder::write_data_structure(const std::string &data_structure_filename) const {
    std::ofstream data_structure_out(data_structure_filename.c_str());
    if(!data_structure_out) {
      std::cerr << "Unable to open \"" << data_structure_filename << "\"" << std::endl;
      exit(-1);
    }
    
    write_data_structure(data_structure_out);
  }

  void Forwarder::write_data_structure(std::ofstream &out) const {
    out << "orig_alphabet_size" << std::endl;
    out << orig_alphabet_size << std::endl;
    
    out << "orig_seq_length" << std::endl;
    out << orig_seq_length << std::endl;
    
    out << "nStates2alphabet_size" << std::endl;
    for(std::map<size_t, size_t>::const_iterator it = nStates2alphabet_size.begin(); it != nStates2alphabet_size.end(); ++it)
      out << (*it).first << " " << (*it).second <<std::endl;

    out << "symbol2pair" << std::endl;
    for(std::map<unsigned, s_pair>::const_iterator it = symbol2pair.begin(); it != symbol2pair.end(); ++it)
      out << (*it).first << " " << (*it).second.first << " " << (*it).second.second <<std::endl;
  }  

  void Forwarder::read_data_structure_from_directory(const std::string &directory) {
    std::string full_filename = directory + "/data_structure";

    std::ifstream data_structure_in(full_filename.c_str());
    if(!data_structure_in) {
      std::cerr << "Unable to open \"" << full_filename << "\"" << std::endl;
      exit(-1);
    }

    read_data_structure_from_stream(data_structure_in);
    data_structure_in.close();
  }

  void Forwarder::read_data_structure_from_stream(std::ifstream &in) {
    read_token_or_exit(in, "orig_alphabet_size");
    orig_alphabet_size = read_or_exit<size_t>(in, "original alphabet size");

    read_token_or_exit(in, "orig_seq_length");
    orig_seq_length = read_or_exit<size_t>(in, "original alphabet size");

    read_token_or_exit(in, "nStates2alphabet_size");
    while(!read_token_or_tell(in, "symbol2pair")) {
      size_t no_states = read_or_exit<size_t>(in, "no_states");
      nStates2alphabet_size[no_states] = read_or_exit<size_t>(in, "alphabet_size");
    }
    
    for(size_t i = orig_alphabet_size; i < nStates2alphabet_size.begin()->second; ++i) {
      unsigned pair_symbol  = read_or_exit<unsigned>(in, "pair symbol");
      unsigned left_symbol  = read_or_exit<unsigned>(in, "left symbol");
      unsigned right_symbol = read_or_exit<unsigned>(in, "right symbol");
      symbol2pair[pair_symbol] = s_pair(left_symbol, right_symbol);
    }
  }

  void Forwarder::read_seqs_from_directory(const std::string &directory) {
    std::vector<std::string> filenames;

    DIR *dpdf;
    struct dirent *epdf;
    
    dpdf = opendir(directory.c_str());
    if (dpdf != NULL) {
      while ((epdf = readdir(dpdf))) {
	if(std::strcmp(epdf->d_name, ".") && std::strcmp(epdf->d_name, "..") && std::strcmp(epdf->d_name, ".svn"))
	  filenames.push_back(epdf->d_name);
      }
    } 
    
    for(std::vector<std::string>::iterator it = filenames.begin(); it != filenames.end(); ++it) {
      std::string filename = (*it);
      size_t no_states = (size_t) std::atol(filename.substr(0, filename.length()-4).c_str());
      read_seq_from_directory(directory, no_states);
    }

  }

  void Forwarder::read_seq_from_directory(const std::string &directory, size_t no_states) {
    std::stringstream ss;
    ss << directory << "/" << no_states << ".seq";
    std::string filename = ss.str();
    std::ifstream in(filename.c_str());
    
    if(!in) {
      std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
      exit(-1);
    }
    
    read_seq_from_directory(in, no_states);
    
    in.close();
  }

  void Forwarder::read_seq_from_directory(std::ifstream &in, size_t no_states) {
    std::vector<unsigned> res;
    unsigned tmp = 0;
    while(in >> tmp)
      res.push_back(tmp);
    nStates2seq[no_states] = res;
  }


  void Forwarder::compute_symbol2scale_and_symbol2matrix(Matrix *symbol2matrix, double *symbol2scale, const Matrix &A, const Matrix &B, const size_t alphabet_size) const{
    // compute C matrices for each symbol in the original alphabet
    make_em_trans_probs_array(symbol2scale, symbol2matrix, A, B);
    
    // compute C matrices for each symbol in the extended alphabet
    for(size_t i = orig_alphabet_size; i < alphabet_size; ++i) {
      const s_pair symbol_pair = symbol2pair.find(unsigned(i))->second;
      const unsigned left_symbol = symbol_pair.second; // the multiplication is done in the reverse direction of the sequence
      const unsigned right_symbol  = symbol_pair.first;
      Matrix &left_matrix  = symbol2matrix[left_symbol];
      Matrix &right_matrix = symbol2matrix[right_symbol];
      
      Matrix::blas_mult(left_matrix, right_matrix, symbol2matrix[i]);
      symbol2scale[i] = std::log( symbol2matrix[i].normalize() ) + symbol2scale[left_symbol] + symbol2scale[right_symbol];
    }
  }

  void Forwarder::pthread_compute_symbol2scale_and_symbol2matrix(double *symbol2scale, Matrix *symbol2matrix, 
							     const Matrix &A, const Matrix &B, 
							     const size_t alphabet_size, 
							     const std::vector<ProcessingDevice*> &devices) const{
    // compute C matrices and scales for each symbol in the original alphabet
    make_em_trans_probs_array(symbol2scale, symbol2matrix, A, B);

    // create control object 
    Stage1JobControl control(orig_alphabet_size, alphabet_size, &symbol2pair, devices.size());

    // do computations for the remaining symbols in parallel
    for(unsigned i = 0; i < devices.size(); ++i)
      devices[i]->computeSymbol2ScaleAndSymbol2Matrix(control);

    // join
    for(unsigned i = 0; i < devices.size(); ++i)
      devices[i]->join();
  }

  size_t Forwarder::get_seq_length(size_t no_states) const {
    for(std::map<size_t, std::vector<unsigned> >::const_iterator it = nStates2seq.begin(); it != nStates2seq.end(); it++) {
      if(it->first >= no_states) {
	return (it->second).size();
      }
    }
    return (nStates2seq.rbegin()->second).size();
  }
  
  size_t Forwarder::get_alphabet_size(size_t no_states) const {
    for(std::map<size_t, std::vector<unsigned> >::const_iterator it = nStates2seq.begin(); it != nStates2seq.end(); ++it) {
      if(it->first >= no_states)
	return nStates2alphabet_size.find(it->first)->second;
    }
    return nStates2alphabet_size.rbegin()->second;
  }

  s_pair Forwarder::get_pair(unsigned symbol) const {
    if(symbol2pair.count(symbol) != 0) 
      return symbol2pair.find(symbol)->second;
    return std::pair<unsigned, unsigned>((unsigned)-1, (unsigned)-1);
  }
} // namespace
