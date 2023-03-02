    swig -lua -c++ -Iinclude WaveFile.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WaveFile.so WaveFile_wrap.cxx -lstdc++ -lm -lluajit    
    