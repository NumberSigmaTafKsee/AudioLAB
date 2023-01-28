    swig -lua -c++ -Iinclude RackFXCompBand.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXCompBand.so RackFXCompBand_wrap.cxx -lstdc++ -lm -lluajit    
    