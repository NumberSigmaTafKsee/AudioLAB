swig -Iinclude/Json -lua -c++ src/json.i
gcc -fmax-errors=1 -Iinclude/Json -O2 -fPIC -shared -o json.so src/json_wrap.cxx -lstdc++ /usr/local/lib/libjsoncpp.a -lm -lluajit -ljsoncpp
