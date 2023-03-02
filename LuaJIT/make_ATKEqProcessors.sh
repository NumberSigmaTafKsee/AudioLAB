    swig -lua -c++ -Iinclude ATKEqProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKEqProcessors.so ATKEqProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    