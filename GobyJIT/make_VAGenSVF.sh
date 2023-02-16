    swig -lua -c++ -Iinclude VAGenSVF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAGenSVF.so VAGenSVF_wrap.cxx -lstdc++ -lm -lluajit    
    