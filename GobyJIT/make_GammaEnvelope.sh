    swig -lua -c++ -Iinclude GammaEnvelope.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o GammaEnvelope.so GammaEnvelope_wrap.cxx -lstdc++ -lm -lluajit    
    