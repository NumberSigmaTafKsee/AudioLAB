swig -Iinclude/Std -lua -c++ -I/usr/local/include jsoncpp.i
gcc -Iinclude/Std -Iinclude -O2 -fPIC -shared -o jsoncpp.so jsoncpp_wrap.cxx -lstdc++ -lm -lluajit -L/usr/local/lib -ljsoncpp
