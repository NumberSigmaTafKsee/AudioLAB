    swig -lua -c++ -Iinclude ATKDynamicProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKDynamicProcessors.so ATKDynamicProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    