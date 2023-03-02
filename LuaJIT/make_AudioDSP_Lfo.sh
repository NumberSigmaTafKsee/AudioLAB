    swig -lua -c++ -Iinclude AudioDSP_Lfo.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioDSP_Lfo.so AudioDSP_Lfo_wrap.cxx -lstdc++ -lm -lluajit    
    