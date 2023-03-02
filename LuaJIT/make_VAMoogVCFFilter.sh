    swig -lua -c++ -Iinclude VAMoogVCFFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogVCFFilter.so VAMoogVCFFilter_wrap.cxx -lstdc++ -lm -lluajit    
    