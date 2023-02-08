swig -lua -c++ -Iinclude analog_krajeski_moog.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o analog_krajeski_moog.so \
analog_krajeski_moog_wrap.cxx -lstdc++ -lm -lluajit
