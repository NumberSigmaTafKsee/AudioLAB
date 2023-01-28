swig -lua -c++ -I/usr/local/include/csound src/csound.i
gcc -fmax-errors=1 -I/usr/local/include/csound -O2 -fPIC -march=native -mavx2 -shared -o csound.so src/csound_wrap.cxx -lstdc++ -lm -lluajit -lcsnd6 -lcsound64
