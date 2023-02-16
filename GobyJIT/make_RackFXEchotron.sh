    swig -lua -c++ -Iinclude RackFXEchotron.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXEchotron.so RackFXEchotron_wrap.cxx -lstdc++ -lm -lluajit    
    