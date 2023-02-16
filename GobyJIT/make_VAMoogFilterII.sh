    swig -lua -c++ -Iinclude VAMoogFilterII.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogFilterII.so VAMoogFilterII_wrap.cxx -lstdc++ -lm -lluajit    
    