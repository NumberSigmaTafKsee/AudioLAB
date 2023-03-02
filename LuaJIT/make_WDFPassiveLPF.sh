    swig -lua -c++ -Iinclude WDFPassiveLPF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WDFPassiveLPF.so WDFPassiveLPF_wrap.cxx -lstdc++ -lm -lluajit    
    