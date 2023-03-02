    swig -lua -c++ -Iinclude AmplifierFold.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AmplifierFold.so AmplifierFold_wrap.cxx -lstdc++ -lm -lluajit    
    