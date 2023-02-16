    swig -lua -c++ -Iinclude RackFXReverb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXReverb.so RackFXReverb_wrap.cxx -lstdc++ -lm -lluajit    
    