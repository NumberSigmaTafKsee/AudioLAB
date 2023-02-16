swig -lua -c++ -Iinclude audio_lfo.i
gcc -fmax-errors=1 -std=c++17 -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o audio_lfo.so audio_lfo_wrap.cxx \
-lstdc++ -lm -lluajit
