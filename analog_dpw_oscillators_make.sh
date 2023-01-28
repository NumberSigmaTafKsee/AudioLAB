swig -lua -c++ -Iinclude analog_dpw_oscillators.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I. -O2 -fPIC -march=native -mavx2 -shared -o analog_dpw_oscillators.so \
analog_dpw_oscillators_wrap.cxx -lstdc++ -lm -lluajit
