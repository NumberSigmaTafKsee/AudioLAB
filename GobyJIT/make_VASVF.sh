    swig -lua -c++ -Iinclude VASVF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VASVF.so VASVF_wrap.cxx -lstdc++ -lm -lluajit    
    