swig -lua -c++ -I/usr/local/include/kissfft src/kissfft_real.i
gcc -I/usr/local/include/kissfft -L/usr/local/lib -O2 -march=native -mavx2 -fPIC -shared -o kissfft_real.so src/kissfft_real_wrap.cxx -lstdc++ -lm -lluajit -lkissfft-float
