swig -lua -c++ -Iinclude analog_blits_oscillators.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o analog_blits_oscillators.so \
analog_blits_oscillators_wrap.cxx -lstdc++ -lm -lluajit
