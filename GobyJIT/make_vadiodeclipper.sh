swig -lua -c++ vadiodeclipper.i
gcc -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o vadiodeclipper.so vadiodeclipper_wrap.cxx -lstdc++ -lm -lluajit
