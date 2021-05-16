
* Build procedure
* Getting started
* Using the C++ library
* Using the Python library
* Encoding HMMs
	* HMM example
	* C++ example
	* Python example
* Encoding sequences
* Executables
	* calibrate
	* build_forwarder
	* forward
	* generate_hmm
	* generate_seq
* Contact

# Build procedure: 

To build and install the library, unzip the directory and execute the
following commands in a terminal:

0. 
If you're on Ubuntu, please follow the instructions here to install and build ATLAS: https://gist.github.com/sangheestyle/ca8ef7796aefadad8773
	
You will run into some missing packages along the way, just keep googling and `apt-get`ing until you're through.  Keep all default options if there are any.

1. 

```bash
	$ cd <path to library>/zipHMM-1.0.1/
	$ cmake .
```
2. 

If you're on Ubuntu or the following make gives an error about threading-related calls not being defined to link with, then do this first:

```bash
	$ grep -rl lpthread ./ | xargs sed -i 's/lpthread/pthread/g' 
```
This simply replaces all occrences of the flag "-lpthread" with "-pthread" and the build goes through after that!  It previously did not like the `l` in `-lpthread`.

3. 

```bash
	$ make
	$ bin/calibrate
	$ make test
	$ make install
```

A successful build (the first `make` command above) on Ubuntu server (running on a VirtualBox VM on Windows 10) with the slim desktop UI, should look like this:
![VMBox instance of Unbuntu server with slim desktop and successful build of zipHMMlib.](https://user-images.githubusercontent.com/1606391/118394241-56572280-b5f8-11eb-9e9c-81683e6551b8.png)

At the end, when you run `make install` you'll receive a print out of important locations (for useage of the library & headers), so save it to a file like so:

```
$ make install > important_install_directories.txt
```

Then you can refer back to it without having to run anything.

---

To build in OS X, the Accellerate framework is required (see
https://developer.apple.com/performance/accelerateframework.html). This
is included in the developer tools installed with XCode (see
https://developer.apple.com/xcode/)


---
To build in Linux CMake must be able to find an installation of a BLAS
implementation. For now the CMake script is set up to use Atlas and to
look for it at /com/extra/ATLAS/3.9.84. This will most likely not work
on your machine. You may therefore have to change line 11 in
zipHMM/CmakeLists.txt:

```bash
  set(ATLAS_ROOT "/com/extra/ATLAS/3.9.84")
```

If you are using a different implementation of BLAS than Atlas you
will have to do a few extra simple changes in zipHMM/CMakelists.txt -
look at line 12, 13, 32, 56, 76, 91, 106, 121, 159 and 192.

bin/calibrate finds the optimal number of threads to use in the
parallelized algorithm and saves the number in a file (default is
~/.ziphmm.devices).


# Getting started

Have a look at zipHMM/cpp_example.cpp and zipHMM/python_example.cpp
and try running the following commands from the root directory.

```bash
$ bin/cpp_example
	
$ cd zipHMM/
$ python python_example.py
```

# Using the C++ library

The main class in the library is Forwarder (forwarder.hpp and
forwarder.cpp). Objects of this class represents an observed sequence
that have been preprocessed such that the likelihood of the sequence
can be obtained from a given HMM (pi, A, B) very fast. To build a new
forwarder object just call the empty constructor:

```c
Forwarder();
```

and to read in a sequence call one of the two read_seq methods:

```c
void Forwarder::read_seq(const std::string &seq_filename, const size_t alphabet_size, 
                         std::vector<size_t> nStatesSave, const size_t min_no_eval = 1);
void Forwarder::read_seq(const std::string &seq_filename, const size_t alphabet_size, 
     			 const size_t no_states, const size_t min_no_eval);
void Forwarder::read_seq(const std::string &seq_filename, const size_t alphabet_size, 
                         const size_t min_no_eval = 1)
```

Here seq_filename is the filename of the file containing the observed
sequence, alphabet_size is the size of the alphabet used in the
observed sequence, nStatesSave is a vector indicating the sizes of the
HMMs the user intends to run the Forwarder object on, and min_no_eval
is a guess of the number of times, the preprocessing will be reused
(if unsure about this then leave it out and use the default value). If
nStatesSave contains the vector (2, 4, 8), datastructures obtaining
the fastest evaluation of the forward algorithm for each of the HMM
state space sizes 2, 4, and 8 will be build. If nStatesSave is left
empty a single datastructure obtaining a very fast evaluation of the
forward algorithm for all HMM state space sizes will be saved.

The second constructor serves as a convenient way to call the first
constructor with only one HMM size in nStatesSave.
The third constructor serves as a convenient way to call the first 
constructor with an empty nStates2save vector.

After building an Forwarder object, it can be saved to disk using the method

```c
void write_to_directory(const std::string &directory) const;
```

Here directory should contain the path (relative to the root of the
library) of the intended location of the datastructure.

To read a previously saved datastructure, one of the following two methods 
can be used:

```c
void Forwarder::read_from_directory(const std::string &dirname);
void Forwarder::read_from_directory(const std::string &directory, const size_t no_states);
```

Using the first one, the entire datastructure is being rebuilt. Using
the second one only the datastructure matching no_states is being
rebuild. This will be faster in many cases. If you did not save the
datastructure for the size of your HMM, then use the first
constructor. The forward algorithm will figure out which of the saved
data structures is most optimal for your HMM.

Finally, to get the loglikelihood of the observed sequence in a
specific model, one of the following methods are used:

```c
double Forwarder::forward(const Matrix &pi, const Matrix &A, const Matrix &B) const;
double Forwarder::pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B, 
                                  const std::string &device_filename = DEFAULT_DEVICE_FILENAME) const;
```

The second method is a parallelized version of the forward algorithm,
whereas the first one is single-threaded. pi, A and B specifies the
HMM parameters. They can either be read from a file or build in C++ as
described below in the section 'Encoding HMMs'. The parallelized version
takes an additional filename as parameter. This filename should be the
path to the file created by the calibrate program, which finds the
optimal number of threads to use in the parallelized forward
algorithm. The default filename is ~/.ziphmm.devices. If you did not
move the file, then leave the device_filename parameter out.

See zipHMM/cpp_example.cpp for a simple example.

# Using the Python library

To use the Python library in another project, copy zipHMM/pyZipHMM.py
and zipHMM/libpyZipHMM.so to the root of your project folder after
building the library and import pyZipHMM. See zipHMM/python_example.py
and zipHMM/python_test.py for details on how to use the library.

A Forwarder object can be constructed from an observed sequence in the
following ways:

```python
from pyZipHMM import *
f = Forwarder.fromSequence(seqFilename = "example.seq", alphabetSize = 3, minNoEvals = 500)
```

To save the datastructure to disk do as follows:

```python
f.writeToDirectory("example_preprocessed")
```

To read a previously saved datastructure from disk use either of the
two methods:

```python
f2 = Forwarder.fromDirectory(directory = "../example_out")
f2 = Forwarder.fromDirectory(directory = "../example_out", nStates = 3)
```

Finally, to evaluate the loglikelihood of the sequence in a given model
(matrices pi, A and B) use either of

```python
loglikelihood = f.forward(pi, A, B)
loglikelihood = f.ptforward(pi, A, B)
```

where the second method is parallelized. The three matrices pi, A and B
can be read from a file or build in Python as described below.

See zipHMM/python_example.py for an example.

# Encoding HMMs

An HMM consists of three matrices: 

   - pi, containing initial state probabilities: pi_i is the
     probability of the model starting in state i.

   - A, containing transition probabilities: A_{ij} is the
   probability of the transition from state i to state j.

   - B, containing emission probabilities: B_{io} is the probability
     of state i emitting symbol o.

These three matrices can either be build in the code (in C++ or
Python) or they can be encoded in a text file. The format of the text
file is as follows:

# HMM example 
---
```
no_states
3
alphabet_size
4
pi
0.1
0.2
0.7
A
0.1 0.2 0.7
0.3 0.4 0.3
0.5 0.5 0.0
B
0.1 0.2 0.3 0.4
0.2 0.3 0.4 0.1
0.3 0.4 0.1 0.2
```
---

To read and write HMMs from and to files in C++, use the methods

```c
void read_HMM(Matrix &resultInitProbs, Matrix &resultTransProbs, Matrix &resultEmProbs, const std::string &filename);
void write_HMM(const Matrix &initProbs, const Matrix &transProbs, const Matrix &emProbs, const std::string &filename);
```

To read and write HMMs from and to files in Python, use the functions

```python
readHMM(filename) -> (pi, A, B)
writeHMM(pi, A, B, filename) -> None
```

To build a matrix in C++ do as illustrated in the following example:

## C++ example

```c
#include "matrix.hpp"

size_t nRows = 3;
size_t nCols = 4;
zipHMM::Matrix m(3,4);
m(0,0) = 0.1;
m(0,1) = 0.2;
...
m(2,3) = 0.2;
```

To build a matrix in Python do as illustrated here:

## Python example

```python
import pyZipHMM

nRows = 3
nCols = 4
m = pyZipHMM.Matrix(nRows, nCols)
m[0,0] = 0.1
m[0,1] = 0.2
...
m[2,3] = 0.2
```

# Encoding sequences

The alphabet of observables are encoded using integers. Thus if the
size of the alphabet is M, the observables are encoded using 0, 1, 2,
..., M - 1. A sequence of observables is encoded in a text file with
the single observations are seperated by whitespace. See example.seq
for an example.


# Executables

## calibrate

Usage: `bin/calibrate`

Finds the optimal number of threads to use in the parallelized version
of the forward algorithm.

## build_forwarder

Usage: `bin/build_forwarder -s <sequence filename> -M <alphabet size> -o <output directory> [-N <number of states>]*`

Builds a Forwarder object from the sequence in the file specified in
<sequence filename> and writes it to the directory specified in
<output directory>. <alphabet size> should be the size of the alphabet
used in the observed sequence, and the file specified in <sequence
filename> should contain a single line containing white space
separated integers between 0 and <alphabet size> - 1. The list of HMM
sizes to generate the data structure for can be specified using the -N parameter.

Examples:
```bash
bin/build_forwarder -s example.seq -M 3 -o example_out
bin/build_forwarder -s example.seq -M 3 -o example_out -N 2
bin/build_forwarder -s example.seq -M 3 -o example_out -N 2 -N 4 -N 8 -N 16
```

## forward

Usage: `bin/forward (-s <sequence filename> -m <HMM filename> [-e #expected forward calls] [-o <output directory>] ) | (-d <preprocessing directory> -m <HMM filename>) [-p]`

Runs the forward algorithm and outputs the loglikelhood. This
executable can be called in two different ways:

```bash
bin/forward -s example.seq -m example.hmm -e 500 -o example_out
bin/forward -d example_out/ -m example.hmm
```

In the first example the loglikelihood is evaluated based on the
observed sequence in example.seq and the HMM specified in
example.hmm. In the second example the loglikelihood is evaluated
based on the previously saved data structure in example_out/ and the
HMM specified in example.hmm. In both cases the -p parameter can be
used to use the parallelized version. In the first example the user
can optionally choose to save the data structure in eg. example_out/
using the -o parameter:

```bash
bin/forward -s example.seq -m example.hmm -e 500 -o example_out/
```

## generate_hmm

Usage: `bin/generate_hmm <number of states> <alphabet size> <HMM filename>`

Generates a random HMM with <number of states> states and <alphabet
size> observables, and saves it to <HMM filename>.

## generate_seq

Usage: `bin/generate_seq <HMM filename> <length> <observed sequence output filename> <state sequence output filename>`

Given an HMM specified in <HMM filename>, runs the HMM for <length>
iterations and saves the resulting sequence of observables to
<observed sequence output filename> and the resulting sequence of
hidden states to <state sequence output filename>.

# Contact

If you encounter any problems or have questions about using this
software, please post an issue here:
	https://github.com/enjoysmath/ziphmm/issues	
	  
	  
# Origination / Inventor Contact Info

This repository is a direct fork of the inventors' code repository: https://github.com/mailund/ziphmm

On that page is contact information for them.

