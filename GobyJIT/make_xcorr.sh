swig -Iinclude -lua -c++ src/xcorr.i
gcc -Iinclude -O2 -march=native -mavx2 -fPIC -shared -o xcorr.so src/xcorr_wrap.cxx src/xcorr.c -lstdc++ -lm -lluajit -lfftw3
