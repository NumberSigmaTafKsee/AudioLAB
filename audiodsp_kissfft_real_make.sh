swig -lua -c++ -I/usr/local/include/kissfft Kits/AudioDSP/audiodsp_kissfft_real.i
gcc -I/usr/local/include/kissfft -L/usr/local/lib -O2 -march=native -mavx2 -fPIC -shared \
-o kissfft_real.so Kits/AudioDSP/audiodsp_kissfft_real_wrap.cxx -lstdc++ -lm -lluajit -lkissfft-float
