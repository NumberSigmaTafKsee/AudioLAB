    swig -lua -c++ -Iinclude audiodsp_envelope.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_envelope.so audiodsp_envelope_wrap.cxx -lstdc++ -lm -lluajit 
    