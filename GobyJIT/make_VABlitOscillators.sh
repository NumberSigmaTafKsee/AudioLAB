    swig -lua -c++ -Iinclude VABlitOscillators.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VABlitOscillators.so VABlitOscillators_wrap.cxx -lstdc++ -lm -lluajit    
    