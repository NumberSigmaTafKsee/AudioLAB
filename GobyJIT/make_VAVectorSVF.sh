    swig -lua -c++ -Iinclude VAVectorSVF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAVectorSVF.so VAVectorSVF_wrap.cxx -lstdc++ -lm -lluajit    
    