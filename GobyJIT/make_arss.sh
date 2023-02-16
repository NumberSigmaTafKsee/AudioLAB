swig -lua -c++ -Isrc src/ARSS.i
gcc -IARSS/src -O2 -fPIC -march=native -mavx2 -shared -o arss.so src/ARSS_wrap.cxx -lstdc++ -lm -lluajit -lfftw3
