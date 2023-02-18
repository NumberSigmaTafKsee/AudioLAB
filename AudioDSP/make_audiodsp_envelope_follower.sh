    swig -lua -c++ -Iinclude audiodsp_envelope_follower.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_envelope_follower.so audiodsp_envelope_follower_wrap.cxx -lstdc++ -lm -lluajit 
    