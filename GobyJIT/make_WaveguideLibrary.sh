    swig -lua -c++ -Iinclude WaveguideLibrary.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WaveguideLibrary.so WaveguideLibrary_wrap.cxx -lstdc++ -lm -lluajit    
    