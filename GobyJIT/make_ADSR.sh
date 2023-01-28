    swig -lua -c++ -Iinclude ADSR.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ADSR.so ADSR_wrap.cxx -lstdc++ -lm -lluajit    
    