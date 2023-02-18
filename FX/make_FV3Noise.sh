    swig -lua -c++ -I../include -I../include FV3Noise.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3Noise.so FV3Noise_wrap.cxx -lstdc++ -lm -lluajit 
    