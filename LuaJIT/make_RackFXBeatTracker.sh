    swig -lua -c++ -Iinclude RackFXBeatTracker.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXBeatTracker.so RackFXBeatTracker_wrap.cxx -lstdc++ -lm -lluajit    
    