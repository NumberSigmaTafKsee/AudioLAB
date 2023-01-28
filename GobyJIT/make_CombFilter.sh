    swig -lua -c++ -Iinclude CombFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o CombFilter.so CombFilter_wrap.cxx -lstdc++ -lm -lluajit    
    