    swig -lua -c++ -Iinclude ToneStack.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ToneStack.so ToneStack_wrap.cxx -lstdc++ -lm -lluajit    
    