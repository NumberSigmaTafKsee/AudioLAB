    swig -lua -c++ -Iinclude VAKorg35LPFFilter.i
    gcc -fmax-errors=1 -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAKorg35LPFFilter.so VAKorg35LPFFilter_wrap.cxx -lstdc++ -lm -lluajit    
    
