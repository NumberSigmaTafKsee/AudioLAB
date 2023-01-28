    swig -lua -c++ -Iinclude VADiodeClipper.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VADiodeClipper.so VADiodeClipper_wrap.cxx -lstdc++ -lm -lluajit    
    