swig -lua -c++ -I/usr/local/include/kissfft Kits/AudioDSP/audiodsp_kissfft.i
gcc -I/usr/local/include/kissfft -L/usr/local/lib -O2 -march=native -mavx2 -fPIC \
-shared -o kissfft.so Kits/AudioDSP/audiodsp_kissfft_wrap.cxx -lstdc++ -lm -lluajit -lkissfft-float
