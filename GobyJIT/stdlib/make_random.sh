swig -Iinclude -lua -c++ random.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -shared -o random.so random_wrap.cxx -lstdc++ -lm -lluajit
