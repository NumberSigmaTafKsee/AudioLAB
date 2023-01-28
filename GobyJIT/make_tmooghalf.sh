swig -lua -c++ -Iinclude src/TMoogHalfLadder.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o TMoogHalfLadder.so src/TMoogHalfLadder_wrap.cxx -lstdc++ -lm -lluajit
