    swig -lua -c++ -Iinclude StereoVector.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o StereoVector.so StereoVector_wrap.cxx -lstdc++ -lm -lluajit    
    