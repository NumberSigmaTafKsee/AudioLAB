swig -Iinclude -lua -c++ scalar.i
gcc -Iinclude -O2 -fPIC -shared -o scalar.so scalar_wrap.cxx -lstdc++ -lm -lluajit
