    swig -lua -c++ -Iinclude VAVecSVF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAVecSVF.so VAVecSVF_wrap.cxx -lstdc++ -lm -lluajit    
    