swig -Iinclude/netsockets -lua -c++ netsockets.i
gcc -Iinclude/netsockets -O2 -fPIC -shared -o netsockets.so netsockets_wrap.cxx lib/liblib_netsockets.a  -lstdc++ -lm -lluajit -lrt -lpthread
