    swig -lua -c++ -Iinclude ClipSerpentCurve.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ClipSerpentCurve.so ClipSerpentCurve_wrap.cxx -lstdc++ -lm -lluajit    
    