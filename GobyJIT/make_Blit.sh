    swig -lua -c++ -Iinclude Blit.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Blit.so Blit_wrap.cxx -lstdc++ -lm -lluajit    
    