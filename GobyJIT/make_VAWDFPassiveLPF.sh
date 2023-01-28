    swig -lua -c++ -Iinclude VAWDFPassiveLPF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAWDFPassiveLPF.so VAWDFPassiveLPF_wrap.cxx -lstdc++ -lm -lluajit    
    