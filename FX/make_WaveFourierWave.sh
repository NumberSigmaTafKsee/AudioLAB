    swig -lua -c++ -I../include -I../include WaveFourierWave.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o WaveFourierWave.so WaveFourierWave_wrap.cxx -lstdc++ -lm -lluajit 
    