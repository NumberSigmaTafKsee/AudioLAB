    swig -lua -c++ -Iinclude StereoDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o StereoDelay.so StereoDelay_wrap.cxx -lstdc++ -lm -lluajit    
    