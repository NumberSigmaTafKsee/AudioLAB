    swig -lua -c++ -Iinclude DistortionProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DistortionProcessors.so DistortionProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    