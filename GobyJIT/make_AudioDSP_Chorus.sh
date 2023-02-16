    swig -lua -c++ -Iinclude AudioDSP_Chorus.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioDSP_Chorus.so AudioDSP_Chorus_wrap.cxx -lstdc++ -lm -lluajit    
    