swig -lua -c++ -Iinclude src/TOberheim.i
gcc -fmax-errors=1 -Iinclude -O2 -march=native -mavx2 -fPIC -shared -o TOberheim.so src/TOberheim_wrap.cxx -lstdc++ -lm -lluajit
