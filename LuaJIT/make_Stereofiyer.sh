    swig -lua -c++ -Iinclude Stereofiyer.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Stereofiyer.so Stereofiyer_wrap.cxx -lstdc++ -lm -lluajit    
    