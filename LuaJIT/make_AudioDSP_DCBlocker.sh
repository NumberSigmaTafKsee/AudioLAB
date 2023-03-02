    swig -lua -c++ -Iinclude AudioDSP_DCBlocker.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioDSP_DCBlocker.so AudioDSP_DCBlocker_wrap.cxx -lstdc++ -lm -lluajit    
    