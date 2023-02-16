    swig -lua -c++ -Iinclude AudioDSP_Delay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioDSP_Delay.so AudioDSP_Delay_wrap.cxx -lstdc++ -lm -lluajit    
    