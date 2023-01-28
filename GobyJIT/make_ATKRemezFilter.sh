    swig -lua -c++ -Iinclude ATKRemezFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKRemezFilter.so ATKRemezFilter_wrap.cxx -lstdc++ -lm -lluajit    
    