swig -lua -c++ src/map.i
gcc -fmax-errors=1  -O2 -fPIC -shared -o map.so src/map_wrap.cxx -lstdc++ -lm -lluajit
