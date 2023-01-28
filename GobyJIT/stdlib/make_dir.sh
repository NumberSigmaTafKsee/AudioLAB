swig -lua -c++ dir.i
gcc -O2 -fPIC -shared -o dir.so dir_wrap.cxx -lstdc++ -lm -lluajit
