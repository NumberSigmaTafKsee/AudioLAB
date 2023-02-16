    swig -lua -c++ -Iinclude OrfanidisEQ.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o OrfanidisEQ.so OrfanidisEQ_wrap.cxx -lstdc++ -lm -lluajit    
    