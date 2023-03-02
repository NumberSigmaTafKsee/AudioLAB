    swig -lua -c++ -Iinclude DynamicProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DynamicProcessors.so DynamicProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    