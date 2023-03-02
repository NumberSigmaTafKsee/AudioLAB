    swig -lua -c++ -Iinclude VCAProcessor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VCAProcessor.so VCAProcessor_wrap.cxx -lstdc++ -lm -lluajit    
    