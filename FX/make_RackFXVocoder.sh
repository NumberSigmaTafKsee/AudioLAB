    swig -lua -c++ -I../include -I../include RackFXVocoder.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXVocoder.so RackFXVocoder_wrap.cxx -lstdc++ -lm -lluajit 
    