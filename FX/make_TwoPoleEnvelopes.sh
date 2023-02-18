    swig -lua -c++ -I../include -I../include TwoPoleEnvelopes.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o TwoPoleEnvelopes.so TwoPoleEnvelopes_wrap.cxx -lstdc++ -lm -lluajit 
    