    swig -lua -c++ -I../include -I../include FXRoboVerb.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FXRoboVerb.so FXRoboVerb_wrap.cxx -lstdc++ -lm -lluajit 
    