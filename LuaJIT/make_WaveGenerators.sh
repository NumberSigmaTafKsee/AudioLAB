    swig -lua -c++ -Iinclude WaveGenerators.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WaveGenerators.so WaveGenerators_wrap.cxx -lstdc++ -lm -lluajit    
    