    swig -lua -c++ -Iinclude SchroederAllpass.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o SchroederAllpass.so SchroederAllpass_wrap.cxx -lstdc++ -lm -lluajit    
    