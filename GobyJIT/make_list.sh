    swig -lua -c++ -Iinclude list.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o list.so list_wrap.cxx -lstdc++ -lm -lluajit    
    