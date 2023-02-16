swig -lua -c++ -Iinclude src/chorus.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o chorus.so src/chorus_wrap.cxx -lstdc++ -lm -lluajit
