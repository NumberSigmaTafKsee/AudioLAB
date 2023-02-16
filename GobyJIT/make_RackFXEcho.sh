    swig -lua -c++ -Iinclude RackFXEcho.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXEcho.so RackFXEcho_wrap.cxx -lstdc++ -lm -lluajit    
    