swig -lua -c++ -Iinclude src/TImprovedMoog.i
gcc -Iinclude -fmax-errors=1 -Iinclude/Filters -O2 -march=native -mavx2 -fPIC -shared -o TImprovedMoog.so src/TImprovedMoog_wrap.cxx -lstdc++ -lm -lluajit
