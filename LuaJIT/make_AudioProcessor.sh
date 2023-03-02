    swig -lua -c++ -Iinclude AudioProcessor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioProcessor.so AudioProcessor_wrap.cxx -lstdc++ -lm -lluajit    
    