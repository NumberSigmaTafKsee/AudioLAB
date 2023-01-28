swig -lua -c++ src/aafilter.i
gcc -I../include -Isrc -O2 -fPIC -march=native -mavx2 -shared -o aafilter.so src/aafilter_wrap.cxx -lstdc++ -lm -lluajit
