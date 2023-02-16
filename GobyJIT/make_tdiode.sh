swig -lua -c++ -Iinclude src/TDiodeLadder.i
gcc -Iinclude -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o TDiodeLadder.so src/TDiodeLadder_wrap.cxx -lstdc++ -lm -lluajit
