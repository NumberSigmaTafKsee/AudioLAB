swig -lua -c++ -Iinclude src/MoogLadders.i
gcc -Iinclude -O2 -fPIC -mavx2 -march=native -shared -o MoogLadders.so src/MoogLadders_wrap.cxx -lstdc++ -lm -lluajit
