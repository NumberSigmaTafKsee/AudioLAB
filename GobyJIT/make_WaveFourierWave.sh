    swig -lua -c++ -Iinclude WaveFourierWave.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WaveFourierWave.so WaveFourierWave_wrap.cxx -lstdc++ -lm -lluajit    
    