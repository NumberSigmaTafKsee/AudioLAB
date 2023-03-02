    swig -lua -c++ -Iinclude BezierDistortion.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o BezierDistortion.so BezierDistortion_wrap.cxx -lstdc++ -lm -lluajit    
    