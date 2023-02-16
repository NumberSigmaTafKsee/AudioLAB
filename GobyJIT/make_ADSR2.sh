    swig -lua -c++ -Iinclude ADSR2.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ADSR2.so ADSR2_wrap.cxx -lstdc++ -lm -lluajit    
    