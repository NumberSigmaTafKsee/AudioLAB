swig -lua -c++ -Iinclude sstwaveshaper.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o sstwaveshaper.so sstwaveshaper_wrap.cxx  -lstdc++ -lm -lluajit
