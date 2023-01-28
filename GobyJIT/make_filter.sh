swig -lua -c++ -Iinclude src/Filter.i
gcc -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Filter.so src/Filter_wrap.cxx -lstdc++ -lm -lluajit
