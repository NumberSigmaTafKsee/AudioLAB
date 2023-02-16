    swig -lua -c++ -Iinclude AudioDSP_FirstOrderAllPass.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioDSP_FirstOrderAllPass.so AudioDSP_FirstOrderAllPass_wrap.cxx -lstdc++ -lm -lluajit    
    