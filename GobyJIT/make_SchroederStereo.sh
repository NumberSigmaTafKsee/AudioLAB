    swig -lua -c++ -Iinclude SchroederStereo.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o SchroederStereo.so SchroederStereo_wrap.cxx -lstdc++ -lm -lluajit    
    