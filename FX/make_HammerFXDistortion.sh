    swig -lua -c++ -I../include -I../include HammerFXDistortion.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o HammerFXDistortion.so HammerFXDistortion_wrap.cxx -lstdc++ -lm -lluajit 
    