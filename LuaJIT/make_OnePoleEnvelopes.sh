    swig -lua -c++ -Iinclude OnePoleEnvelopes.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o OnePoleEnvelopes.so OnePoleEnvelopes_wrap.cxx -lstdc++ -lm -lluajit    
    