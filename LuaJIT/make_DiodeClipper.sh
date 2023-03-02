    swig -lua -c++ -Iinclude DiodeClipper.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DiodeClipper.so DiodeClipper_wrap.cxx -lstdc++ -lm -lluajit    
    