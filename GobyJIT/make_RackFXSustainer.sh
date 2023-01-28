    swig -lua -c++ -Iinclude RackFXSustainer.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXSustainer.so RackFXSustainer_wrap.cxx -lstdc++ -lm -lluajit    
    