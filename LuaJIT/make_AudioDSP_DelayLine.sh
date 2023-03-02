    swig -lua -c++ -Iinclude AudioDSP_DelayLine.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioDSP_DelayLine.so AudioDSP_DelayLine_wrap.cxx -lstdc++ -lm -lluajit    
    