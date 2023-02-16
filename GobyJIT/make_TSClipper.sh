    swig -lua -c++ -Iinclude TSClipper.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o TSClipper.so TSClipper_wrap.cxx -lstdc++ -lm -lluajit    
    