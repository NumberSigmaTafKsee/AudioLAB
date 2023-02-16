swig -lua -c++ -Iinclude analog_microtracker_filter.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I. -O2 -fPIC -march=native -mavx2 -shared -o analog_microtracker_filter.so \
analog_microtracker_filter_wrap.cxx -lstdc++ -lm -lluajit
