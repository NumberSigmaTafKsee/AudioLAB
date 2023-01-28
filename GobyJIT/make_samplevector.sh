swig -lua -c++ -Iinclude SampleVector.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o sv.so SampleVector_wrap.cxx -lstdc++ -lm -lluajit
