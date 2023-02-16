    swig -lua -c++ -Iinclude RackFXConvolotron.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXConvolotron.so RackFXConvolotron_wrap.cxx -lstdc++ -lm -lluajit    
    