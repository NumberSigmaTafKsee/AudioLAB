    swig -lua -c++ -Iinclude ATKBlockLMSFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKBlockLMSFilter.so ATKBlockLMSFilter_wrap.cxx -lstdc++ -lm -lluajit    
    