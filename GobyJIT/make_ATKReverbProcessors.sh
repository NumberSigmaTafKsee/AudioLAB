    swig -lua -c++ -Iinclude ATKReverbProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKReverbProcessors.so ATKReverbProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    