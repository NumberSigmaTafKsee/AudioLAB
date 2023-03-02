    swig -lua -c++ -Iinclude DiodeSim.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DiodeSim.so DiodeSim_wrap.cxx -lstdc++ -lm -lluajit    
    