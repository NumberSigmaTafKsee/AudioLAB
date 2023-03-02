    swig -lua -c++ -Iinclude DCFProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DCFProcessors.so DCFProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    