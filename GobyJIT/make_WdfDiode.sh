    swig -lua -c++ -Iinclude WdfDiode.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WdfDiode.so WdfDiode_wrap.cxx -lstdc++ -lm -lluajit    
    