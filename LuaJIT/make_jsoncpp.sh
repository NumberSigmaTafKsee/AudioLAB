swig -Iinclude/Std -lua -c++ -I/usr/local/include src/jsoncpp.i
gcc -Iinclude/Std -Iinclude -O2 -fPIC -shared -o jsoncpp.so src/jsoncpp_wrap.cxx -lstdc++ -lm -lluajit -L/usr/local/lib -ljsoncpp
