swig -lua -c++ map.i
gcc -fmax-errors=1  -O2 -fPIC -shared -o map.so map_wrap.cxx -lstdc++ -lm -lluajit
