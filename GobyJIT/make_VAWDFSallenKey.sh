    swig -lua -c++ -Iinclude VAWDFSallenKey.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAWDFSallenKey.so VAWDFSallenKey_wrap.cxx -lstdc++ -lm -lluajit    
    