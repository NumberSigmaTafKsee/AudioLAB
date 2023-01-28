    swig -lua -c++ -Iinclude RandomSmoother.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RandomSmoother.so RandomSmoother_wrap.cxx -lstdc++ -lm -lluajit    
    