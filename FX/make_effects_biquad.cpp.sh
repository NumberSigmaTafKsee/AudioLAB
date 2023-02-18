    swig -lua -c++ -I../include -I../include effects_biquad.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o effects_biquad.cpp.so effects_biquad.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    