    swig -lua -c++ -Iinclude RackFXVocoder.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXVocoder.so RackFXVocoder_wrap.cxx -lstdc++ -lm -lluajit    
    