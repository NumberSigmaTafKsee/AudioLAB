    swig -lua -c++ -Iinclude ChebyDistortion.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ChebyDistortion.so ChebyDistortion_wrap.cxx -lstdc++ -lm -lluajit    
    