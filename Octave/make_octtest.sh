swig -octave -c++ octave_test.i
mkoctfile -o octave_test octave_test_wrap.cxx -lstdc++ -lm
