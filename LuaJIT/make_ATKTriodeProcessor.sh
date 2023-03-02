    swig -lua -c++ -Iinclude ATKTriodeProcessor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKTriodeProcessor.so ATKTriodeProcessor_wrap.cxx -lstdc++ -lm -lluajit    
    