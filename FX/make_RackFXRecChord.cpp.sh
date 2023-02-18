    swig -lua -c++ -I../include -I../include RackFXRecChord.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXRecChord.cpp.so RackFXRecChord.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    