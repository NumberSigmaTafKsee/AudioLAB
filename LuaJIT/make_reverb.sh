    swig -lua -c++ -Iinclude reverb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o reverb.so reverb_wrap.cxx -lstdc++ -lm -lluajit    
    