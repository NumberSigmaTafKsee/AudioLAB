    swig -lua -c++ -Iinclude FxDSPMIDI.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPMIDI.so FxDSPMIDI_wrap.cxx -lstdc++ -lm -lluajit    
    