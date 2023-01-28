    swig -lua -c++ -Iinclude KHSchreoderReverb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o KHSchreoderReverb.so KHSchreoderReverb_wrap.cxx -lstdc++ -lm -lluajit    
    