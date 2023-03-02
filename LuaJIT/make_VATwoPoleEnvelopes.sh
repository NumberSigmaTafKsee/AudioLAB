    swig -lua -c++ -Iinclude VATwoPoleEnvelopes.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VATwoPoleEnvelopes.so VATwoPoleEnvelopes_wrap.cxx -lstdc++ -lm -lluajit    
    