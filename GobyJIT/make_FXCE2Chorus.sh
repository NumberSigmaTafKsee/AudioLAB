    swig -lua -c++ -Iinclude FXCE2Chorus.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXCE2Chorus.so FXCE2Chorus_wrap.cxx -lstdc++ -lm -lluajit    
    