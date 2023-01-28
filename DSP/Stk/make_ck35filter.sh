swig -lua -c++ -Iinclude ck35filter.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o ck35filter.so ck35filter_wrap.cxx  -lstdc++ -lm -lluajit
