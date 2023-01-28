    swig -lua -c++ -Iinclude PolygonalOscillator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o PolygonalOscillator.so PolygonalOscillator_wrap.cxx -lstdc++ -lm -lluajit    
    