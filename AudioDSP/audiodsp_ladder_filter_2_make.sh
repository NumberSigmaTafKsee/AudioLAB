swig -lua -c++ -Iinclude analog_ladder_filter_2.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I. -O2 -fPIC -march=native -mavx2 -shared -o analog_ladder_filter_2.so \
analog_ladder_filter_2_wrap.cxx -lstdc++ -lm -lluajit
