    swig -lua -c++ -Iinclude WDFSallenKey.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WDFSallenKey.so WDFSallenKey_wrap.cxx -lstdc++ -lm -lluajit    
    