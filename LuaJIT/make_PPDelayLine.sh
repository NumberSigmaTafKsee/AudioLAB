    swig -lua -c++ -Iinclude PPDelayLine.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o PPDelayLine.so PPDelayLine_wrap.cxx -lstdc++ -lm -lluajit    
    