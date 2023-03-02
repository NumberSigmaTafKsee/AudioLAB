    swig -lua -c++ -Iinclude ATKLMSFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKLMSFilter.so ATKLMSFilter_wrap.cxx -lstdc++ -lm -lluajit    
    