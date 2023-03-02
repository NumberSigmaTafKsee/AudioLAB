swig -lua -c++ -Iinclude src/TWaveTable.i
gcc -Iinclude -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o TWaveTable.so src/TWaveTable_wrap.cxx -lstdc++ -lm -lluajit
