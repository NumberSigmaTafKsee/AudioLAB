swig -lua -c++ src/kissfft.i
gcc -I/usr/local/include/kissfft -L/usr/local/lib -O2 -march=native -mavx2 -fPIC -shared -o kissfft.so src/kissfft_wrap.cxx -lstdc++ -lm -lluajit -lkissfft-float
