swig -Iinclude -lua -c++ src/random.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -O2 -fPIC -shared -o random.so src/random_wrap.cxx -lstdc++ -lm -lluajit
