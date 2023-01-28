swig -lua -c++ -Iinclude cppfilters.i
gcc -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o cppfilters.so cppfilters_wrap.cxx -lstdc++ -lm -lluajit
