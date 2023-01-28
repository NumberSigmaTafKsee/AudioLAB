swig -lua -c++ -Iinclude src/TMoogLadder.i
gcc -fmax-errors=1 -Iinclude -O2 -march=native -mavx2 -fPIC -shared -o TMoogLadder.so src/TMoogLadder_wrap.cxx -lstdc++ -lm -lluajit
