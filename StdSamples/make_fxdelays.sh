swig -lua -c++ -Iinclude audiodsp_fxdelays.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o fxdelays.so audiodsp_fxdelays_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
