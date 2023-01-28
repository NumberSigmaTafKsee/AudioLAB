    swig -lua -c++ -Iinclude VAMS20Filter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMS20Filter.so VAMS20Filter_wrap.cxx -lstdc++ -lm -lluajit    
    