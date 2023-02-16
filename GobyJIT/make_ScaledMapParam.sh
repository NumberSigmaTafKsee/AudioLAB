    swig -lua -c++ -Iinclude ScaledMapParam.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ScaledMapParam.so ScaledMapParam_wrap.cxx -lstdc++ -lm -lluajit    
    