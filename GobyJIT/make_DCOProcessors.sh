    swig -lua -c++ -Iinclude DCOProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DCOProcessors.so DCOProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    