swig -Iinclude -lua -c++ hashmap.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -shared -o hashmap.so  hashmap_wrap.cxx -lstdc++ -lm -lluajit
