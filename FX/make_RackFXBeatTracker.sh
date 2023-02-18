    swig -lua -c++ -I../include -I../include RackFXBeatTracker.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXBeatTracker.so RackFXBeatTracker_wrap.cxx -lstdc++ -lm -lluajit 
    