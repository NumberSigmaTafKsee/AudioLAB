    swig -lua -c++ -Iinclude VAPolygonalOscillator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAPolygonalOscillator.so VAPolygonalOscillator_wrap.cxx -lstdc++ -lm -lluajit    
    