    swig -lua -c++ -Iinclude DSF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DSF.so DSF_wrap.cxx -lstdc++ -lm -lluajit    
    