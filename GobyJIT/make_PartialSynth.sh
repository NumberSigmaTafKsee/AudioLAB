    swig -lua -c++ -Iinclude PartialSynth.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o PartialSynth.so PartialSynth_wrap.cxx -lstdc++ -lm -lluajit    
    