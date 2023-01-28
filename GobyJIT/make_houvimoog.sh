swig -lua -c++ -Iinclude src/THouvilainenMoog.i
gcc -Iinclude -fmax-errors=1 -Iinclude/Filters -O2 -march=native -mavx2 -fPIC -shared -o THouvilainenMoog.so src/THouvilainenMoog_wrap.cxx -lstdc++ -lm -lluajit
