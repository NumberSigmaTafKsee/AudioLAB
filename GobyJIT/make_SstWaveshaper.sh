    swig -lua -c++ -Iinclude SstWaveshaper.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o SstWaveshaper.so SstWaveshaper_wrap.cxx -lstdc++ -lm -lluajit    
    