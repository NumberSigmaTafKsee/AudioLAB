    swig -lua -c++ -Iinclude ATKDistortionProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKDistortionProcessors.so ATKDistortionProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    