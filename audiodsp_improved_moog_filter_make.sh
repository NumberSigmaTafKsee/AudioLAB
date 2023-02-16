swig -lua -c++ -Iinclude analog_improved_moog_filter.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I. -O2 -fPIC -march=native -mavx2 -shared -o analog_improved_moog_filter.so \
analog_improved_moog_filter_wrap.cxx -lstdc++ -lm -lluajit
