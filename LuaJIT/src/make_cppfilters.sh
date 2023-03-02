swig -lua -c++ -Iinclude src/cppfilters.i
gcc -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o cppfilters.so src/cppfilters_wrap.cxx -lstdc++ -lm -lluajit
