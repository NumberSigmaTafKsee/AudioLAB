swig -lua -c++ -Iinclude cdca.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o cdca.so cdca_wrap.cxx  -lstdc++ -lm -lluajit
