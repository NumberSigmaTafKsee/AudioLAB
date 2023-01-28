    swig -lua -c++ -Iinclude BasicDelayLineStereo.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o BasicDelayLineStereo.so BasicDelayLineStereo_wrap.cxx -lstdc++ -lm -lluajit    
    