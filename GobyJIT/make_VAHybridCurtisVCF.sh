    swig -lua -c++ -Iinclude VAHybridCurtisVCF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAHybridCurtisVCF.so VAHybridCurtisVCF_wrap.cxx -lstdc++ -lm -lluajit    
    