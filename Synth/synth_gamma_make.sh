swig -lua -c++ -Iinclude synth_gamma.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o gamma.so synth_gamma_wrap.cxx lib/libGamma.a -lstdc++ -lm -lportaudio -lluajit
