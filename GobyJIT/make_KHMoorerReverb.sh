    swig -lua -c++ -Iinclude KHMoorerReverb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o KHMoorerReverb.so KHMoorerReverb_wrap.cxx -lstdc++ -lm -lluajit    
    