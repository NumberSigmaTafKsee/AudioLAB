    swig -lua -c++ -I../include -I../include RackFXOpticalTrem.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXOpticalTrem.so RackFXOpticalTrem_wrap.cxx -lstdc++ -lm -lluajit 
    