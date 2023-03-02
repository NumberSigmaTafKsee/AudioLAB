    swig -lua -c++ -Iinclude CEnvelopeDetector.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o CEnvelopeDetector.so CEnvelopeDetector_wrap.cxx -lstdc++ -lm -lluajit    
    