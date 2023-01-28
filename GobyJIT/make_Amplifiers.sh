    swig -lua -c++ -Iinclude Amplifiers.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Amplifiers.so Amplifiers_wrap.cxx -lstdc++ -lm -lluajit    
    