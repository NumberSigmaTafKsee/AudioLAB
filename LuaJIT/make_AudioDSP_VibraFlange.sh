    swig -lua -c++ -Iinclude AudioDSP_VibraFlange.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioDSP_VibraFlange.so AudioDSP_VibraFlange_wrap.cxx -lstdc++ -lm -lluajit    
    