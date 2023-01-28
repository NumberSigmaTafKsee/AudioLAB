swig -lua -c++ vablitoscillator.i
gcc -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o vablitoscillator.so vablitoscillator_wrap.cxx -lstdc++ -lm -lluajit
