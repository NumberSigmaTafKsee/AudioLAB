    swig -lua -c++ -I../include -I../include RackFXMusicDelay.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXMusicDelay.so RackFXMusicDelay_wrap.cxx -lstdc++ -lm -lluajit 
    