    swig -lua -c++ -Iinclude TwoPoleEnvelopes.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o TwoPoleEnvelopes.so TwoPoleEnvelopes_wrap.cxx -lstdc++ -lm -lluajit    
    