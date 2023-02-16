swig -lua -c++ -I/usr/local/include/libnyquist nyquist.i
gcc -I/usr/local/include/libnyquist -O2 -fPIC -shared -o nyquist.so nyquist_wrap.cxx -lstdc++ -lm -lluajit -llibnyquist
