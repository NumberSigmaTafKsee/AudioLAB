    swig -lua -c++ -I../include -I../include AudioDSP_VibraFlange.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioDSP_VibraFlange.so AudioDSP_VibraFlange_wrap.cxx -lstdc++ -lm -lluajit 
    