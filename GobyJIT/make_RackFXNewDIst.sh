    swig -lua -c++ -Iinclude RackFXNewDIst.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXNewDIst.so RackFXNewDIst_wrap.cxx -lstdc++ -lm -lluajit    
    