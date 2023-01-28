    swig -lua -c++ -Iinclude Old.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Old.so Old_wrap.cxx -lstdc++ -lm -lluajit    
    