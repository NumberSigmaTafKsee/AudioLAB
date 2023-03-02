    swig -lua -c++ -Iinclude FXRoboVerb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXRoboVerb.so FXRoboVerb_wrap.cxx -lstdc++ -lm -lluajit    
    