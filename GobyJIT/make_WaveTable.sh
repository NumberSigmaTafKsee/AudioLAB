    swig -lua -c++ -Iinclude WaveTable.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WaveTable.so WaveTable_wrap.cxx -lstdc++ -lm -lluajit    
    