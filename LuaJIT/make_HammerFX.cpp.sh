    swig -lua -c++ -Iinclude HammerFX.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFX.cpp.so HammerFX.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    