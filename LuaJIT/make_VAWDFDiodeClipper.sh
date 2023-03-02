    swig -lua -c++ -Iinclude VAWDFDiodeClipper.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAWDFDiodeClipper.so VAWDFDiodeClipper_wrap.cxx -lstdc++ -lm -lluajit    
    