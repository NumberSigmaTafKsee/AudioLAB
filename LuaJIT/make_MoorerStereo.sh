    swig -lua -c++ -Iinclude MoorerStereo.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o MoorerStereo.so MoorerStereo_wrap.cxx -lstdc++ -lm -lluajit    
    