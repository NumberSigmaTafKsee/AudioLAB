    swig -lua -c++ -Iinclude PeakHoldCascade.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o PeakHoldCascade.so PeakHoldCascade_wrap.cxx -lstdc++ -lm -lluajit    
    