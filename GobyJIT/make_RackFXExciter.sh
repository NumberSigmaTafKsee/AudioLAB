    swig -lua -c++ -Iinclude RackFXExciter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXExciter.so RackFXExciter_wrap.cxx -lstdc++ -lm -lluajit    
    