    swig -lua -c++ -Iinclude OnePoleEnvelope.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o OnePoleEnvelope.so OnePoleEnvelope_wrap.cxx -lstdc++ -lm -lluajit    
    