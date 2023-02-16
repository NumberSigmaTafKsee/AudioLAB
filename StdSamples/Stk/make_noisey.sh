swig -Iinclude -lua -c++ src/noisey.i
gcc -Iinclude -fmax-errors=1 -O2 -fPIC -shared -o noisey.so src/noisey_wrap.cxx -lstdc++ -lm -lluajit
