swig -lua -c++ -Iinclude audiodsp_fxobjects.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o fxobjects.so audiodsp_fxobjects_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
