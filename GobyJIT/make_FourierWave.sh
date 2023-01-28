    swig -lua -c++ -Iinclude FourierWave.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FourierWave.so FourierWave_wrap.cxx -lstdc++ -lm -lluajit    
    