    swig -lua -c++ -Iinclude RackFXOpticalTrem.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXOpticalTrem.so RackFXOpticalTrem_wrap.cxx -lstdc++ -lm -lluajit    
    