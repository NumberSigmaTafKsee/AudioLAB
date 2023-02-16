    swig -lua -c++ -Iinclude FXProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXProcessors.so FXProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    