swig -lua -c++ src/fir_filter.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o fir_filter.so src/fir_filter_wrap.cxx src/fir_filter.cpp -lstdc++ -lm -lluajit
