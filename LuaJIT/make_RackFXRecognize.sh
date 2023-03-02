    swig -lua -c++ -Iinclude RackFXRecognize.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXRecognize.so RackFXRecognize_wrap.cxx -lstdc++ -lm -lluajit    
    