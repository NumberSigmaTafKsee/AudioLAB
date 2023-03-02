swig -lua -c++ cpumemory.i
gcc -fmax-errors=1 -O2 -march=native -mavx2 -fPIC -shared -o cpumemory.so cpumemory_wrap.cxx cpumemory.cpp -lstdc++ -lm -lluajit
